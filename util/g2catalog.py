#
# g2catalog.py -- Fits viewer interface to the Gen2 star catalog
# 
# Eric Jeschke (eric@naoj.org)
#
import remoteObjects as ro

from ginga.util import wcs
from ginga.util import catalog as StarCatalog
from ginga.misc import Bunch, Task

from Gen2.starlist import starlist
from Gen2.starlist import starfilter


class CatalogServer(object):

    def __init__(self, logger, full_name, key, dbhost, description):
        ## super(CatalogServer, self).__init__(logger, full_name, key, None,
        ##                                         description)
        self.logger = logger
        self.full_name = full_name
        self.short_name = key
        self.description = description
        self.kind = 'g2starcatalog'
        # {name} {ra} {dec} {mag} {flag} {b_r} {preference} {priority} {dst}
        self.index = { 'name': 'name',
                       'ra': 'ra', 'dec': 'dec',
                       'mag': 'mag', 'flag': 'flag',
                       'field': 'field',
                       'b-r': 'b_r', 'preference': 'preference',
                       'priority': 'priority', 'dst': 'dst',
                       }
        self.format = 'deg'
        self.equinox = 2000.0
        self.dbhost = dbhost
        self.threadPool = None
        self.svcname = dbhost

        # For compatibility with URL catalog servers
        self.params = {}
        count = 0
        for label, key in (('RA', 'ra'), ('DEC', 'dec'),
                           ('Min Radius', 'r1'), ('Max Radius', 'r2'),
                           ('Max Mag', 'm1'), ('Min Mag', 'm2'),
                           ('Catalog', 'catalog')):
            self.params[key] = Bunch.Bunch(name=key, convert=str,
                                           label=label, order=count)
            count += 1

        #self.reset_conn()

    def getParams(self):
        return self.params

    def set_index(self, **kwdargs):
        self.index.update(kwdargs)

    def reset_conn(self, threadPool):
        self.threadPool = threadPool

        # add anything here for one time setup
        self.catalog = starlist.CatalogSearch(dbhost=self.dbhost, logger=self.logger, threadpool=self.threadPool)

    def search(self, **params):
        """For compatibility with 'regular' star catalogs."""

        print "search params=%s" % params
        ra, dec = params['ra'], params['dec']
        if not (':' in ra):
            # Assume RA and DEC are in degrees
            ra_deg = float(ra)
            dec_deg = float(dec)
        else:
            # Assume RA and DEC are in standard string notation
            ra_deg = wcs.hmsStrToDeg(ra)
            dec_deg = wcs.dmsStrToDeg(dec)

        # Subaru uses degrees for fov and no minimum radius
        fov_deg = float(params['r2']) / 60.0

        catalog = params['catalog'].strip()
        if not catalog:
            catalog = "usnob,gsc,sao"

        # Default min and max magnitudes if none specified
        s = params['m1'].strip()
        if not s:
            lowermag = 0.0
        else:
            lowermag = float(s)

        s = params['m2'].strip()
        if not s:
            uppermag = 21.0
        else:
            uppermag = float(s)
        
        kwdargs = dict(ra=ra_deg, dec=dec_deg, fov=fov_deg, 
                       lowermag=lowermag, uppermag=uppermag,
                       catalog=catalog)
        self.logger.debug("search params are %s" % (str(kwdargs)))


        query_result = self.catalog.search_starcatalog(ra=kwdargs['ra'], dec=kwdargs['dec'], fov=kwdargs['fov'], lowermag=kwdargs['lowermag'], uppermag=kwdargs['uppermag'], catalog=kwdargs['catalog'])
        #print "QUERY RESULT=", query_result
        
        starlist = self.process_starlist(query_result)
        #print "STARLIST=", starlist

        # metadata about the list
        columns = [('Name', 'name'),
                   ('RA', 'ra'),
                   ('DEC', 'dec'),
                   ('Mag', 'mag'),
                   ('Priority', 'priority'),
                   ('Preference', 'preference'),
                   ('Flag', 'flag'),
                   ('DB', 'field'),
                   #('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=1)
        return starlist, info


    def search_sh(self, ra_deg=None, dec_deg=None, equinox=None,
                  fov_deg=None, upper_mag=13.0):
        # For SH:
        # note: limitmag=13.0 is fixed value for sh
        params = dict(ra=ra_deg, dec=dec_deg, fov=fov_deg, equinox=equinox,
                      limitmag=upper_mag)

        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(params)))

        starlist = self.catalog.search_starcatalog(ra=params['ra'], dec=params['dec'], fov=params['fov'], lowermag=0, uppermag=params['limitmag'])

        
        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_sh_stars(kargs=params, starlist=starlist)

        result = { 'selected_stars': starlist, 'prefered_num': len(starlist) }
        return result


    def search_ag(self, ra_deg=None, dec_deg=None, fov_deg=None,
                  probe_ra_deg=None, probe_dec_deg=None, equinox=None,
                  focus=None, inst_name=None,  probe_r=None, probe_theta=None,
                  probe_x=None, probe_y=None, pos_ang_deg=None, upper_mag=None,
                  pref_mag=None, fov_pattern=None):

        # For AG:
        # focus in { CS, NS_IR, NS_OPT, P_OPT, P_IR }
        # inst in { MOIRCS, etc. }
        # note: ra, dec, fov, probe_ra, probe_dec, pa are in degrees
        params = dict(ra=ra_deg, dec=dec_deg, fov=fov_deg, equinox=equinox,
                      probe_ra=probe_ra_deg, probe_dec=probe_dec_deg, 
                      focus=focus, ins=inst_name, probe_r=probe_r, 
                      probe_theta=probe_theta, probe_x=probe_x,
                      probe_y=probe_y, pa=pos_ang_deg, limitmag=upper_mag,
                      goodmag=pref_mag, fov_pattern=fov_pattern)

        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(params)))

        # CHANGE vvvvvv
        starlist = self.catalog.search_starcatalog(ra=params['ra'], dec=params['dec'], fov=params['fov'], lowermag=0, uppermag=params['limitmag'], pa=params['pa'], focus=params['focus'])
     
        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_ag_stars(kargs=params, starlist=starlist)

        result = { 'selected_stars': starlist, 'prefered_num': len(starlist) }         
        return result

    def process_starlist(self, starlist):

        desirable_flags = set((2, 3))
        results = []
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
                if cmp(self.equinox, 2000.0) != 0:
                    ra_deg, dec_deg = wcs.eqToEq2000(ra_deg, dec_deg,
                                                     self.equinox)
                
                ra_txt = wcs.raDegToString(ra_deg, format='%02d:%02d:%06.3f')
                dec_txt = wcs.decDegToString(dec_deg,
                                               format='%s%02d:%02d:%05.2f')
                args['ra'] = ra_txt
                args['dec'] = dec_txt
                args['ra_deg'] = ra_deg
                args['dec_deg'] = dec_deg

                flags = args['flag']
                args['pick'] = (flags in desirable_flags)

                results.append(StarCatalog.Star(**args))

            except Exception, e:
                self.logger.error("Error parsing catalog query results: %s" % (
                    str(e)))
                raise e

        return results

    def process_result(self, query_result):
        #print query_result
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
                   ('DB', 'field'),
                   ('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=query_result['prefered_num'])

        results = self.process_starlist(starlist)
        return info, results

#END
