#
# g2catalog.py -- Fits viewer interface to the Gen2 star catalog
#
# Eric Jeschke (eric@naoj.org)
# T. Inagaki
#

from math import radians
import time
import json
import requests

from numpy import power as pow

from ginga.util import wcs
from ginga.util import catalog as StarCatalog
from ginga.misc import Bunch, Task

from Gen2.starlist import starlist
from Gen2.starlist import starfilter

from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u

from astroquery.utils.tap.core import TapPlus
from astroquery.irsa import Irsa
from astroquery.gaia import Gaia as GaiaCatalog
from astropy.coordinates import Angle
from astroquery.vizier import Vizier


TAP_service = TapPlus(url="http://vao.stsci.edu/PS1DR2/tapservice.aspx")


class CatalogServerError(Exception):
    pass


class CatalogServer(object):

    def __init__(self, logger, full_name, key, dbhost, description):
        self.logger = logger
        self.full_name = full_name
        self.short_name = key
        self.description = description
        self.kind = 'g2starcatalog'
        # {name} {ra} {dec} {mag} {flag} {b_r} {preference} {priority} {dst}
        self.index = { 'name': 'name', 'cat_id': 'star_id',
                       'ra': 'ra', 'dec': 'dec',
                       'mag': 'mag', 'flag': 'flag',
                       'catalog': 'field',
                       'b-r': 'b_r', 'preference': 'preference',
                       'priority': 'priority', 'dst': 'dst',
                       'description': 'description',
                       }
        self.format = 'deg'
        self.equinox = 2000.0
        self.dbhost = dbhost
        self.threadPool = None
        self.svcname = dbhost

        # For compatibility with URL catalog servers
        self.params = {}
        count = 0
        for label, key in (('RA', 'ra'), ('DEC', 'dec'), ('Equinox', 'equinox'),
                           ('Min Radius', 'r1'), ('Max Radius', 'r2'),
                           ('Max Mag', 'm1'), ('Min Mag', 'm2'),
                           ('Pos Angle', 'pa'), ('Focus', 'focus'),
                           ('Catalog', 'catalog')):
            self.params[key] = Bunch.Bunch(name=key, convert=str,
                                           label=label, order=count)
            count += 1

    def get_params(self):
        return self.params

    getParams = get_params

    def set_index(self, **kwdargs):
        self.index.update(kwdargs)

    def reset_conn(self, threadPool):
        self.threadPool = threadPool

        # add anything here for one time setup
        self.catalog = starlist.CatalogSearch(dbhost=self.dbhost,
                                              logger=self.logger,
                                              threadpool=self.threadPool)

    def get_search_params(self, params):

        # "stringify" all params to make compatible with the GUI
        # fields
        params = dict(list(map(lambda item: (item[0], str(item[1])),
                            params.items())))

        self.logger.debug('params={}'.format(params))
        ra, dec = params['ra'], params['dec']
        if not (':' in ra):
            # Assume RA and DEC are in degrees
            ra_deg = float(ra)
            dec_deg = float(dec)
        else:
            # Assume RA and DEC are in standard string notation
            ra_deg = wcs.hmsStrToDeg(ra)
            dec_deg = wcs.dmsStrToDeg(dec)

        s = params.get('equinox', '2000.0').strip()
        if not s:
            equinox = 2000.0
        else:
            equinox = float(s)

        # Subaru uses degrees for fov and no minimum radius
        fov_deg = float(params['r2']) / 60.0

        # Get PA
        s = params.get('pa', '').strip()
        if not s:
            pa = 0.0
        else:
            pa = float(s)

        # Get foci
        focus = params.get('focus', '').strip()

        catname = params.get('catalog', '').strip()
        self.logger.debug('catname={}'.format(catname))

        if catname.upper() == 'SUBARU':
            catalog = "usnob,gsc,sao"
        elif catname.upper() == 'PANSTARRS':
            catalog = "panstarrs"
        elif catname.upper() == 'UCAC4':
            catalog = "ucac4_sources"
        elif catname.upper() == 'GAIA_WEB':
            catalog = "gaiadr2.gaia_source"
        # this is for SH and HSC catalog query
        elif catname.upper() == '2MASS':
            catalog = "II/246"
        elif catname ==  '':
            catname = "subaru"
            catalog = "usnob,gsc,sao"
        else:
            #catname = "SUBARU"
            #catalog = "usnob,gsc,sao"
            catalog = '{}'.format(catname)

        self.logger.info('catname={}  catalog={}'.format(catname, catalog))

        # Default min and max magnitudes if none specified
        s = params.get('m1', '').strip()
        if not s:
            lowermag = 0.0
        else:
            lowermag = float(s)

        s = params.get('m2', '').strip()
        if not s:
            uppermag = 21.0
        else:
            uppermag = float(s)

        kwdargs = dict(ra=ra_deg, dec=dec_deg, equinox=equinox,
                       fov=fov_deg, catname=catname,
                       lowermag=lowermag, uppermag=uppermag,
                       catalog=catalog, pa=pa, focus=focus)
        self.logger.debug("search params are %s" % (str(kwdargs)))
        return kwdargs

    def search(self, **params):
        kwdargs = self.get_search_params(params)

        # TODO: temporary assign hard-coded value 0 as lower mag
        starlist = self.catalog.search_starcatalog(kwdargs['ra'], kwdargs['dec'],
                                                   kwdargs['fov'], 0.0, #kwdargs['lowermag'],
                                                   kwdargs['uppermag'],
                                                   catalog=kwdargs['catalog'])

        starlist = self.process_starlist(starlist)

        # metadata about the list
        columns = [('Name', 'name'),
                   ('RA', 'ra'),
                   ('DEC', 'dec'),
                   ('Mag', 'mag'),
                   ('Priority', 'priority'),
                   ('Preference', 'preference'),
                   ('Flag', 'flag'),
                   ('DB', 'catalog'),
                   ('ID', 'cat_id'),
                   ('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=1)
        return starlist, info


    def process_starlist(self, starlist):

        desirable_flags = set((2, 3))
        results = []
        reassigned = []
        self.logger.debug("Iterating over %d results" % (
            len(starlist)))
        for elts in starlist:
            try:
                # Extract particulars of this catalog's format into
                # our known keywords
                args = {}
                for key, idx in self.index.items():
                    #args[key] = elts[idx]
                    args[key] = elts.get(idx, None)

                # Standardize on the convention for RA/DEC.  ra/dec are in
                # traditional notation as strings and ra_deg/dec_deg are
                # floats
                if (self.format == 'deg') or not (':' in args['ra']):
                    # Assume RA and DEC are in degrees
                    ra_deg = float(args['ra'])
                    dec_deg = float(args['dec'])
                else:
                    # Assume RA and DEC are in standard string notation
                    ra_deg = wcs.hmsStrToDeg(args['ra'])
                    dec_deg = wcs.dmsStrToDeg(args['dec'])

                # convert ra/dec via EQUINOX change if catalog EQUINOX is
                # not the same as our default one (2000)
                if int(self.equinox) != 2000:
                    ra_deg, dec_deg = wcs.eqToEq2000(ra_deg, dec_deg,
                                                     self.equinox)

                ra_txt = wcs.ra_deg_to_str(ra_deg)
                dec_txt = wcs.dec_deg_to_str(dec_deg)
                args['ra'] = ra_txt
                args['dec'] = dec_txt
                args['ra_deg'] = ra_deg
                args['dec_deg'] = dec_deg

                flags = args['flag']
                args['pick'] = (flags in desirable_flags)

                star = StarCatalog.Star(**args)

                # adjust priority if star is blocklisted
                if blocklist.check_blocklist(star):
                    star['priority'] = 9999999
                    star['preference'] = 9999999
                    star['description'] = 'BLOCKLISTED'
                    reassigned.append(star)
                else:
                    results.append(star)

            except Exception as e:
                self.logger.error("Error parsing catalog query results: %s" % (
                    str(e)))
                raise e

        results.extend(reassigned)
        return results

    def process_result(self, query_result):
        starlist = query_result['selected_stars']

        # filter stars here
        # TODO

        # metadata about the list
        columns = [('Name', 'name'),
                   ('RA', 'ra'),
                   ('DEC', 'dec'),
                   ('Mag', 'mag'),
                   ('Priority', 'priority'),
                   ('Preference', 'preference'),
                   ('Flag', 'flag'),
                   ('DB', 'catalog'),
                   ('ID', 'cat_id'),
                   ('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=query_result['prefered_num'])

        results = self.process_starlist(starlist)
        return info, results


class WebCatalogs(CatalogServer):

    def web_search(self, search_params):
        self.logger.debug('webcatalog server....')

        catname = search_params['catname']
        self.logger.debug('catname={}'.format(catname))

        if catname.upper() == 'PANSTARRS':
            self.logger.debug("start Panstarrs3.searching...")
            starlist = PanStarrs3.search(search_params, self.logger)
        elif catname.upper() == 'UCAC4':
            self.logger.debug("start Ucac4.searching...")
            starlist = Ucac4.search(search_params, self.logger)
        elif catname.upper() == 'GAIA_WEB':
            self.logger.debug("start Gaia.searching...")
            starlist = Gaia.search(search_params, self.logger)
        elif catname.upper() == '2MASS':
            self.logger.debug("start 2mass searching...")
            starlist = TwoMass.search(search_params, self.logger)
        else:
            raise CatalogServerError("error: invalid catalog name={}".format(catname))

        return starlist


class AgWebCatalog(WebCatalogs):

    def search(self, **params):

        self.logger.debug('ag webcatalog server....')

        search_params = self.get_search_params(params)
        self.logger.debug('search_params={}'.format(search_params))


        starlist = self.web_search(search_params)

        #print('starlist={}'.format(starlist))
        filter_params = dict(ra=search_params['ra'], dec=search_params['dec'],
                             equinox=params['equinox'], fov=search_params['fov'],
                             pa=search_params['pa'],
                             probe_ra=params['probe_ra_deg'],
                             probe_dec=params['probe_dec_deg'],
                             focus=search_params['focus'], ins=params['inst_name'],
                             probe_r=params['probe_r'],
                             probe_theta=params['probe_theta'],
                             probe_x=params['probe_x'],
                             probe_y=params['probe_y'],
                             limitmag=search_params['uppermag'],
                             goodmag=search_params['lowermag'],
                             fov_pattern=params['fov_pattern'])

        #print('filter_params={}'.format(filter_params))
        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_ag_stars(filter_params, starlist)

        #print('star filter={}'.format(starlist))
        query_result = { 'selected_stars': starlist,
                         'query_params': search_params,
                         'prefered_num': len(starlist) }

        return query_result


class ShWebCatalog(WebCatalogs):

    def search(self, **params):

        self.logger.debug('sh webcatalog server....')

        search_params = self.get_search_params(params)
        self.logger.debug('search_params={}'.format(search_params))

        # overwrite lower/upper mag, assign targetmag
        search_params['lowermag'] = params['lowermag']
        search_params['uppermag'] = params['uppermag']
        search_params['targetmag'] = params['targetmag']

        starlist = self.web_search(search_params)

        self.logger.info("catalog search returned {} stars".format(len(starlist)))

        # Filter stars for SH:
        filter_params = dict(ra=search_params['ra'], dec=search_params['dec'],
                             equinox=params['equinox'], fov=search_params['fov'],
                             pa=search_params['pa'],
                             #focus=k['focus'],
                             limitmag=search_params['uppermag'],
                             goodmag=search_params['targetmag']
                             )

        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_sh_stars(filter_params, starlist)

        query_result = { 'selected_stars': starlist,
                         'query_params': search_params,
                         'prefered_num': len(starlist) }
        return query_result


class HscWebCatalog(WebCatalogs):

    def search(self, **params):

        self.logger.debug('hsc webcatalog server....')

        search_params = self.get_search_params(params)
        self.logger.debug('search_params={}'.format(search_params))

        starlist = self.web_search(search_params)

        starlist = self.process_starlist(starlist)
        #print "STARLIST=", starlist

        # metadata about the list
        columns = [('Name', 'name'),
                   ('RA', 'ra'),
                   ('DEC', 'dec'),
                   ('Mag', 'mag'),
                   ('Priority', 'priority'),
                   ('Preference', 'preference'),
                   ('Flag', 'flag'),
                   ('DB', 'catalog'),
                   ('ID', 'cat_id'),
                   ('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=1)
        return starlist, info


class TwoMass:

    @staticmethod
    def search(params, logger):
        logger.info(f'2Mass params={params}')

        starlist = []
        append = starlist.append

        catalog = params['catalog']
        lower_mag = params['lowermag']
        upper_mag = params['uppermag']

        v = Vizier(columns=["RAJ2000", "DEJ2000", "2MASS", "Qflg", "Cflg", "Jmag", "Hmag", "Kmag"], catalog=catalog, row_limit=-1, column_filters={'Qflg': '=AAA', 'Cflg': '=000', 'Kmag': f'{lower_mag}..{upper_mag}'})
        stars = v.query_region(SkyCoord(ra=params['ra'], dec=params['dec'], unit=(u.deg, u.deg)), radius=Angle(params['fov'], "deg"))[0]

        for num in range(len(stars)):
            #logger.info(f"ra={star['RAJ2000']}, dec={star['DEJ2000']}, name={star['_2MASS']}, mag={star['Kmag']}")
            try:
                ra = stars['RAJ2000'][num]
                dec = stars['DEJ2000'][num]
                ra_rad = radians(ra)
                dec_rad = radians(dec)
                star_id = stars['_2MASS'][num]
                name = f"2MASS{star_id}"
                mag = stars['Kmag'][num]
            except Exception as e:
                logger.error(f'key error={e}')
                continue
            else:
                pass
                append(dict(star_id=star_id, name=name, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=mag, flag=2,  b_r=99.9, r_i=99.9, field=catalog))

        logger.info(f'starlist={len(starlist)}')

        return starlist


class Gaia:

    def flag(parallax):
    #def flag(pmra, pmdec, parallax):

        """
           Star Classification.
           Algorithm is designed by Akihito Tajitsu (Subaru Telescope, NAOJ)
           flag 2, the best guiding star candidate
           flag 5, might be ok, but not preferred

            algorithm below is commented out now.
            but might be used later if parallax algorithm is not working well
        """
        # try:
        #     if (pmra+pmdec) > 10000: # find out the number(10000?)
        #         flag = 5
        #         return flag
        # except Exception as e:
        #     print('error. {}'.format(e))
        #     flag = 5
        #     return flag

        try:
            assert parallax > 0.01
            flag = 2
        except Exception as e:
            flag = 5

        return flag


    @staticmethod
    def search(params, logger):
        logger.warning('Gaia params={}'.format(params))

        starlist = []
        append = starlist.append

        lowermag = 0
        params['lowermag'] = lowermag

        query_time = time.time()

        # note: can limit the query result.  e.g.,  select top 500 ...
        sql = """ SELECT designation, source_id, ra, dec, phot_rp_mean_mag, parallax, pmra, pmdec FROM {catalog} WHERE CONTAINS(POINT('ICRS', {catalog}.ra, {catalog}.dec),CIRCLE('ICRS',{ra},{dec},{fov}))=1 and phot_rp_mean_mag > {lowermag} and phot_rp_mean_mag <= {uppermag}; """.format(**params)

        logger.debug('Gaia sql={}'.format(sql))
        job = GaiaCatalog.launch_job_async(sql, dump_to_file=False)
        logger.debug('Gaia query done={}'.format(time.time()-query_time))

        catalog = params['catalog']
        gaias = job.get_results()

        for gaia in gaias:
            ra = gaia['ra']
            dec = gaia['dec']
            ra_rad = radians(ra)
            dec_rad = radians(dec)
            # py37 is ok, but py38 spits decode error
            #name = gaia['designation'].decode('utf-8')
            name = gaia['designation']
            append(dict(star_id=gaia['source_id'], name=name, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=gaia['phot_rp_mean_mag'], flag=Gaia.flag(gaia['parallax']),  b_r=99.9, r_i=99.9, field=catalog))

        return starlist


class PanStarrs3:

    def flag(rMag, rKronMag, iMag, iKronMag):

        """
           star classification. algorithm is designed by Ichi Tanaka(Subaru Telescope)
           flag 2, the best guiding star candidate
           flag 3, good candiate
           flag 5, not preferred
           flag 10, avoid the candidate

        """

        """ note: cosider this algorithm provided by Tajitsu san. talk to Ichi san.
            for now, commented out. use algorithm below for more detail star classification.
        if rApMag - rPSFMag  > -0.1
            flag = 2
        else:
            flag = 5
        """

        lower = -0.5
        upper = 0.2

        rmag_diff =  rMag - rKronMag
        imag_diff = iMag - iKronMag

        #print('rmag_diff={}, imag_diff={}'.format(rmag_diff, imag_diff))

        if not rmag_diff:
            rmag_diff = -1

        if not imag_diff:
            imag_diff = -1

        if (lower <= rmag_diff <= upper) and (lower <= imag_diff <= upper):
            flag = 2
        elif (lower <= rmag_diff <= upper) or (lower <= imag_diff <= upper):
            flag = 3
        elif (rmag_diff > upper) and (imag_diff > upper):
            flag = 5
        elif (rmag_diff < lower) and (imag_diff < lower):
            flag = 10
        else:
            flag = 10

        #print('rmag_diff={}, imag_diff={}, flag={}'.format(rmag_diff, imag_diff, flag))

        return flag

    @staticmethod
    def search(params, logger):
        logger.debug('Panstarrs3 params={}'.format(params))

        lowermag = 0
        #print('im here... params={}'.format(params))
        uppermag = params['uppermag']

        constraints = {"nDetections.gt": 5,
                       "rPSFMag.lte": uppermag, "rPSFMag.gt": lowermag,
                       "iPSFMag.lte": uppermag, "iPSFMag.gt": lowermag}

        data = constraints.copy()

        data['ra'] = params['ra']
        data['dec'] = params['dec']
        data['radius'] = params['fov']
        columns = ['objID', 'objName', 'raStack', 'decStack', 'rPSFMag', 'rKronMag', 'iPSFMag', 'iKronMag']
        data['columns'] = '[{}]'.format(','.join(columns))

        baseurl = 'https://catalogs.mast.stsci.edu/api/v0.1'
        release = 'dr2'
        table = 'stack'
        _format = 'json'
        url = "{0}/{1}/{2}/{3}.{4}".format(baseurl, params['catalog'], release, table, _format)

        logger.debug('Panstarrs3 data={}'.format(data))
        logger.debug('Panstarrs3 url={}'.format(url))

        #for _id, name, ra, dec, rmag, rkmag, imag, ikmag in s.get('data'):
        starlist = []
        append = starlist.append

        q_time = time.time()
        res = requests.get(url, params=data)
        logger.debug('Panstarrs3 query done={}'.format(time.time()-q_time))

        #print('url={}'.format(res.url))
        catalog= params['catalog']
        for oid, oname, ra, dec, rmag, rkmag, imag, ikmag in res.json().get('data'):
            ra_rad = radians(ra)
            dec_rad = radians(dec)
            append(dict(star_id=oid, name=oname, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=rmag, flag=PanStarrs.flag(rmag, rkmag, imag, ikmag),  b_r=99.9, r_i=99.9, field=catalog))

        res.close()
        return starlist


class PanStarrs2:

    def flag(rMag, rKronMag, iMag, iKronMag):

        """
           star classification. algorithm is provided by Ichi Tanaka(Subaru Telescope)
           flag 2, the best guiding star candidate
           flag 3, good candiate
           flag 5, not preferred
           flag 10, avoid the candidate

        """

        lower = -0.5
        upper = 0.2

        rmag_diff =  rMag - rKronMag
        imag_diff = iMag - iKronMag

        #print('rmag_diff={}, imag_diff={}'.format(rmag_diff, imag_diff))

        if not rmag_diff:
            rmag_diff = -1

        if not imag_diff:
            imag_diff = -1

        if (lower <= rmag_diff <= upper) and (lower <= imag_diff <= upper):
            flag = 2
        elif (lower <= rmag_diff <= upper) or (lower <= imag_diff <= upper):
            flag = 3
        elif (rmag_diff > upper) and (imag_diff > upper):
            flag = 5
        elif (rmag_diff < lower) and (imag_diff < lower):
            flag = 10
        else:
            flag = 10

        #print('rmag_diff={}, imag_diff={}, flag={}'.format(rmag_diff, imag_diff, flag))

        return flag

    @staticmethod
    def search(params, logger):
        logger.debug('Panstarrs2 params={}'.format(params))

        lowermag = 0.0
        #print('start panstarrs query...')
        #print('params={}'.format(params))

        params['lowermag'] = lowermag
        #print('params w lowermag={}'.format(params))

        query_time = time.time()

        query = """SELECT objID, objName, raStack, decStack, rPSFMag, rKronMag, iPSFMag, iKronMag from dbo.StackObjectView where CONTAINS(POINT('ICRS', raStack, decStack),CIRCLE('ICRS', {ra}, {dec}, {fov}))=1 and nr > 5 and ni > 5 and rPSFMag > {lowermag} and rPSFMag < {uppermag} and iPSFMag > {lowermag} and iPSFMag < {uppermag} and nDetections> 5""".format(**params)
        logger.debug('Panstarrs2 sql={}'.format(query))

        try:
            job = TAP_service.launch_job_async(query)
        except Exception as e:
            logger.error('error: query. {}'.format(e))
            raise CatalogServerError('error: Panstarrs2 query. {}'.format(e))

        logger.debug('panstarrs2 query done={}'.format(time.time()-query_time))

        starlist = []
        append = starlist.append

        catalog = params['catalog']
        for oid, oname, ra, dec, rmag, rkmag, imag, ikmag in job.get_results():
            ra_rad = radians(ra)
            dec_rad = radians(dec)
            append(dict(star_id=oid, name=oname, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=rmag, flag=PanStarrs.flag(rmag, rkmag, imag, ikmag),  b_r=99.9, r_i=99.9, field=catalog))

        return starlist


class Ucac4():

    def flag(le, xx):

        """
           star classification. algorithm is provided by Ichi Tanaka(Subaru Telescope)
           flag 2, the best guiding star candidate
           flag 10, avoid the candidate
        """


        if le == 0 and xx == 0:
            flag = 2
        else:
            flag = 10

        return flag

        # note: might wanna use dsf to improve star classification later.
        #       dsf(combined double star flag)
        #ucac4['dsf'] == 0

    @staticmethod
    def search(params, logger):
        logger.debug('UCAC4 params={}'.format(params))

        c = SkyCoord(ra=params['ra'], dec=params['dec'], unit=(u.deg, u.deg), frame='icrs')

        starlist = []
        append = starlist.append
        query_time = time.time()
        ucac4s = Irsa.query_region(c, catalog=params["catalog"], spatial="Cone", radius= params['fov']*u.deg)

        logger.debug('UCAC4 query done={}'.format(time.time()-query_time))

        lowermag = 0.0
        mask = (ucac4s['rmag'] <= params['uppermag']) & (ucac4s['rmag'] > lowermag) & (ucac4s['imag'] <= params['uppermag']) & (ucac4s['imag'] > lowermag)

        catalog = params['catalog']
        for ucac4 in ucac4s[mask]:
        #for ucac4 in ucac4s:
            name = 'Ucac4{}'.format(ucac4['ucac4_id']) #ucac4['ucac4_id']
            ra = ucac4['ra']
            dec = ucac4['dec']
            ra_rad = radians(ra)
            dec_rad = radians(dec)

            append(dict(star_id=ucac4['uniqueid'], name=name, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=ucac4['rmag'], flag=Ucac4.flag(ucac4['le'], ucac4['xx']),  b_r=99.9, r_i=99.9, field=catalog))


        #print('ucac4={}'.format(starlist))
        return starlist


class PanStarrs:

    def flag(rMag, rKronMag, iMag, iKronMag):

        """
           Star Classification
           Algorithm is provided by Ichi Tanaka (Subaru Telescope, NAOJ)
           flag 2, best guiding star candidate
           flag 3, good candiate
           flag 5, not preferred
           flag 10, avoid the candidate

        """

        lower = -0.5
        upper = 0.2

        rmag_diff =  rMag - rKronMag
        imag_diff = iMag - iKronMag

        #print('rmag_diff={}, imag_diff={}'.format(rmag_diff, imag_diff))

        if not rmag_diff:
            rmag_diff = -1

        if not imag_diff:
            imag_diff = -1

        if (lower <= rmag_diff <= upper) and (lower <= imag_diff <= upper):
            flag = 2
        elif (lower <= rmag_diff <= upper) or (lower <= imag_diff <= upper):
            flag = 3
        elif (rmag_diff > upper) and (imag_diff > upper):
            flag = 5
        elif (rmag_diff < lower) and (imag_diff < lower):
            flag = 10
        else:
            flag = 10

        #print('rmag_diff={}, imag_diff={}, flag={}'.format(rmag_diff, imag_diff, flag))

        return flag

    @staticmethod
    def search(params, logger):
        logger.debug('Panstarrs params={}'.format(params))

        data_release = 'dr2' # options: dr1 is older, dr2 is the latest and best
        table = 'stacked' # options: 'mean'. for our purpose, stacked is better option

        c = SkyCoord(ra=params['ra'], dec=params['dec'], unit=(u.deg, u.deg), frame='icrs')

        starlist = []
        append = starlist.append

        if table == "stacked":
            print('stacked')

            query_time = time.time()
            cat = Catalogs.query_region(c, params['fov'], catalog=params['catalog'],
                                        data_release=data_release, table=table)
            logger.debug('panstarrs query done. time={}'.format(time.time()-query_time))

            # note: lower mag is fixed value 0
            lowermag = 0.0
            mask = (cat['rPSFMag'] <= params['uppermag']) & (cat['rPSFMag'] > lowermag) & (cat['iPSFMag'] <= params['uppermag']) & (cat['iPSFMag'] > lowermag) & ((cat['nr'] > 5) | (cat['ni'] > 5))

            """ note:
                nr: Number of single epoch detections in r filter.
                ni: Number of single epoch detections in i filter.
            """

            catalog = params['catalog']
            for star in cat[mask]:
                _id = star['objID']
                name = star['objName']
                ra = star['raStack']
                dec = star['decStack']
                ra_rad = radians(ra)
                dec_rad = radians(dec)
                mag = star['rPSFMag']
                append(dict(star_id=_id, name=name, ra_rad=ra_rad, dec_rad=dec_rad, ra=ra, dec=dec, mag=mag, flag=PanStarrs.flag(star['rPSFMag'], star['rKronMag'], star['iPSFMag'], star['iKronMag']),  b_r=99.9, r_i=99.9, field=catalog))


        return starlist


class AgCatalogServer(CatalogServer):

    def search(self, **params):

        k = self.get_search_params(params)

        # Query the catalog
        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(k)))

        # TODO: temporary assign hard-coded value 0 as lower mag
        starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
                                                   k['fov'], 0, #k['lowermag'],
                                                   k['uppermag'],
                                                   pa=k['pa'],
                                                   focus=k['focus'],
                                                   catalog=k['catalog'])

        self.logger.debug("catalog search returned %d stars" % (len(starlist)))

        # Filter stars for AG:
        # focus in { CS, NS_IR, NS_OPT, P_OPT, P_IR }
        # inst in { MOIRCS, etc. }
        # note: ra, dec, fov, probe_ra, probe_dec, pa are in degrees
        filter_params = dict(ra=k['ra'], dec=k['dec'],
                             equinox=params['equinox'], fov=k['fov'],
                             pa=k['pa'],
                             probe_ra=params['probe_ra_deg'],
                             probe_dec=params['probe_dec_deg'],
                             focus=k['focus'], ins=params['inst_name'],
                             probe_r=params['probe_r'],
                             probe_theta=params['probe_theta'],
                             probe_x=params['probe_x'],
                             probe_y=params['probe_y'],
                             limitmag=k['uppermag'],
                             goodmag=k['lowermag'],
                             fov_pattern=params['fov_pattern'])

        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_ag_stars(filter_params, starlist)

        query_result = { 'selected_stars': starlist,
                         'query_params': k,
                         'prefered_num': len(starlist) }
        return query_result


class ShCatalogServer(CatalogServer):

    def search(self, **params):

        k = self.get_search_params(params)

        # overwrite lower/upper mag, assign targetmag
        k['lowermag'] = params['lowermag']
        k['uppermag'] = params['uppermag']
        k['targetmag'] = params['targetmag']

        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(k)))
        starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
                                                   k['fov'], k['lowermag'],
                                                   k['uppermag'],
                                                   catalog=k['catalog'])
        self.logger.debug("catalog search returned %d stars" % (len(starlist)))

        # Filter stars for SH:
        filter_params = dict(ra=k['ra'], dec=k['dec'],
                             equinox=params['equinox'], fov=k['fov'],
                             pa=k['pa'],
                             #focus=k['focus'],
                             limitmag=k['uppermag'],
                             goodmag=k['targetmag']
                             )

        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_sh_stars(filter_params, starlist)

        query_result = { 'selected_stars': starlist,
                         'query_params': k,
                         'prefered_num': len(starlist) }
        return query_result


class GuideStarBlocklist(object):

    def __init__(self, filepath):
        # TODO: read in the blocklist from file
        self.filepath = filepath
        self.blocklist = {}
        try:
            with open(self.filepath, 'r') as in_f:
                buf = in_f.read()

            for line in buf.split('\n'):
                line = line.strip()
                # ignore comments and empty lines
                if line.startswith('#') or len(line) == 0:
                    continue
                if '#' in line:
                    fields = line.split('#')
                    line = fields[0]
                    value = fields[1:]
                else:
                    value = True
                cat, cat_id = line.split(',')
                cat_id = int(cat_id)
                self.blocklist[(cat, cat_id)] = value

        except IOError as e:
            pass

    def mk_key(self, star):
        key = (star['catalog'].lower(), star['cat_id'])
        return key

    def check_blocklist(self, star):
        return self.mk_key(star) in self.blocklist

    def checkpoint_file(self, logger):
        # TODO: find a reasonable file format library
        with open(self.filepath, 'w') as out_f:
            for key, value in self.blocklist.items():
                logger.debug(f'key={key}, value={value}')
                catalog, cat_id = key
                out_f.write(f"{catalog},{cat_id}   # {value}\n")

    def add_blocklist(self, star, logger):
        info = (star['ra'], star['dec'], star['name'])
        logger.debug(f'add blocklist. info={info}')
        self.blocklist[self.mk_key(star)] = info
        logger.debug(f'blocklist={blocklist}')
        self.checkpoint_file(logger)

    def remove_blocklist(self, star, logger):
        del self.blocklist[self.mk_key(star)]
        self.checkpoint_file(logger)


# This implements the global blocklist for guide stars
import os
if 'CONFHOME' in os.environ:
    blocklist_path = os.path.join(os.environ['CONFHOME'], 'guideview',
                                  "blocklist.txt")
else:
    blocklist_path = os.path.join("/tmp", "blocklist.txt")

blocklist = GuideStarBlocklist(blocklist_path)

#END
