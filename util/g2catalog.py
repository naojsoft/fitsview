#
# g2catalog.py -- Fits viewer interface to the Gen2 star catalog
# 
# Eric Jeschke (eric@naoj.org)
#
from ginga.util import wcs
from ginga.util import catalog as StarCatalog
from ginga.misc import Bunch, Task

from Gen2.starlist import starlist
from Gen2.starlist import starfilter


class CatalogServer(object):

    def __init__(self, logger, full_name, key, dbhost, description):
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
        for label, key in (('RA', 'ra'), ('DEC', 'dec'), ('Equinox', 'equinox'),
                           ('Min Radius', 'r1'), ('Max Radius', 'r2'),
                           ('Max Mag', 'm1'), ('Min Mag', 'm2'),
                           ('Pos Angle', 'pa'), ('Focus', 'focus'),
                           ('Catalog', 'catalog')):
            self.params[key] = Bunch.Bunch(name=key, convert=str,
                                           label=label, order=count)
            count += 1

    def getParams(self):
        return self.params

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
        params = dict(map(lambda item: (item[0], str(item[1])),
                          params.items()))
        
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

        catalog = params.get('catalog', '').strip()
        if not catalog:
            catalog = "usnob,gsc,sao"

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
                       fov=fov_deg, 
                       lowermag=lowermag, uppermag=uppermag,
                       catalog=catalog, pa=pa, focus=focus)
        self.logger.debug("search params are %s" % (str(kwdargs)))
        return kwdargs

    def search(self, **params):
        kwdargs = self.get_search_params(params)

        starlist = self.catalog.search_starcatalog(kwdargs['ra'], kwdargs['dec'], kwdargs['fov'], kwdargs['lowermag'], kwdargs['uppermag'], 
                                                   catalog=kwdargs['catalog'])
        #print "QUERY RESULT=", query_result
        
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
                   ('DB', 'field'),
                   #('Description', 'description'),
                   ]
        info = Bunch.Bunch(columns=columns, color='Mag',
                           num_preferred=1)
        return starlist, info


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


class AgCatalogServer(CatalogServer):

    def search(self, **params):

        k = self.get_search_params(params)

        # Query the catalog
        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(k)))
        starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
                                                   k['fov'], k['lowermag'],
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

        self.logger.debug("querying the server %s, params=%s" % (
            self.svcname, str(k)))
        starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
                                                   k['fov'], k['lowermag'],
                                                   k['uppermag'],
                                                   catalog=k['catalog'])
        self.logger.debug("catalog search returned %d stars" % (len(starlist)))
        
        # Filter stars for SH:
        # note: limitmag=13.0 is fixed value for sh
        filter_params = dict(ra=k['ra'], dec=k['dec'],
                             equinox=params['equinox'], fov=k['fov'],
                             pa=k['pa'],
                             #focus=k['focus'], 
                             limitmag=k['uppermag'],
                             #goodmag=k['lowermag']
                             )

        star_select = starfilter.StarSelection(logger=self.logger)
        starlist = star_select.select_sh_stars(filter_params, starlist)

        query_result = { 'selected_stars': starlist,
                         'query_params': k,
                         'prefered_num': len(starlist) }
        return query_result

#END
