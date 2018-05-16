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
from six.moves import map
from six.moves import zip

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

        # starlist = self.catalog.search_starcatalog(kwdargs['ra'], kwdargs['dec'], kwdargs['fov'], kwdargs['lowermag'], kwdargs['uppermag'],
        #                                            catalog=kwdargs['catalog'])

        # TO DO:
        # temporary assign hard-coded value 0 as lower mag
        starlist = self.catalog.search_starcatalog(kwdargs['ra'], kwdargs['dec'], kwdargs['fov'], 0.0, kwdargs['uppermag'],
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

                star = StarCatalog.Star(**args)

                # adjust priority if star is blacklisted
                if blacklist.check_blacklist(star):
                    star['priority'] = 9999999
                    star['preference'] = 9999999
                    star['description'] = 'BLACKLISTED'
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
                   ('DB', 'catalog'),
                   ('ID', 'cat_id'),
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

        # starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
        #                                            k['fov'], k['lowermag'],
        #                                            k['uppermag'],
        #                                            pa=k['pa'],
        #                                            focus=k['focus'],
        #                                            catalog=k['catalog'])

        # TO DO:
        # temporary assign hard-coded value 0 as lower mag
        starlist = self.catalog.search_starcatalog(k['ra'], k['dec'],
                                                   k['fov'], 0,
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


class GuideStarBlacklist(object):

    def __init__(self, filepath):
        # TODO: read in the blacklist from file
        self.filepath = filepath

        self.blacklist = {}
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
                self.blacklist[(cat, cat_id)] = value

        except IOError as e:
            pass

    def mk_key(self, star):
        key = (star['catalog'].lower(), star['cat_id'])
        return key

    def check_blacklist(self, star):
        return self.mk_key(star) in self.blacklist

    def checkpoint_file(self):
        # TODO: find a reasonable file format library
        with open(self.filepath, 'w') as out_f:
            for key, value in self.blacklist.items():
                catalog, cat_id = key
                out_f.write("%s,%d   # %s\n" % (catalog, cat_id, str(value)))

    def add_blacklist(self, star):
        info = (star['ra'], star['dec'], star['name'])
        self.blacklist[self.mk_key(star)] = info
        self.checkpoint_file()

    def remove_blacklist(self, star):
        del self.blacklist[self.mk_key(star)]
        self.checkpoint_file()


# This implements the global blacklist for guide stars
import os
if 'CONFHOME' in os.environ:
    blacklist_path = os.path.join(os.environ['CONFHOME'], 'guideview',
                                  "blacklist.txt")
else:
    blacklist_path = os.path.join("/tmp", "blacklist.txt")

blacklist = GuideStarBlacklist(blacklist_path)

#END
