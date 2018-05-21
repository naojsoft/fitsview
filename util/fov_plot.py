#
# fov_plot.py -- classes implementing graphics to be drawn on the
#                 FOV for guide star selection in Gen2 (see VGW plugin)
#
# Eric Jeschke
# Takeshi Inagaki
#
import math
from operator import itemgetter

from ginga.canvas.CanvasObject import get_canvas_types
from ginga.misc import Bunch
from ginga.util.six.moves import map

cvtypes = get_canvas_types()

import SOSS.GuiderInt.ag_config as ag_config


class TELESCOPEfov(object):
    """A Generic telescope foci FOV object.
    """

    def __init__(self, logger, image, p):
        self.logger = logger
        self.p = p
        ## inner_fov = 0.004167

        # thickness of the marked rings for FOV
        self.ring_thickness = 2
        self.colors = Bunch.Bunch(inst='magenta', outer='red', inner='red',
                                  vignette='green', probe='cyan')

        # calculate radius of probe outer movable area fov
        outer_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, p.outer_fov)
        self.logger.debug("Probe outer movable area radius is %d pixels." % (
            outer_radius))

        ## # calculate radius of probe inner movable area fov
        ## inner_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, inner_fov)
        ## self.logger.debug("Probe inner movable area radius is %d pixels." % (
        ##     inner_radius))

        # calculate probe circle
        probe_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, p.probe_head_fov)
        self.logger.debug("Actual probe head radius is %d pixels." % (
            probe_radius))

        # calculate vignetting mapping fov
        vignette_map = self.calculate_vignette_fov(image, p.ctr_x, p.ctr_y,
                                                   p.f_select, p.ag_pa,
                                                   p.probe_theta, p.probe_r,
                                                   p.fov_pattern)
        self.logger.debug("Vignette map: %s" % (
            str(vignette_map)))

        ## self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.probe_radius = probe_radius
        self.vignette_map = vignette_map


    def calculate_vignette_fov(self, image, ctr_x, ctr_y, focus, pa_deg, theta, r, pattern):
        vlist = ag_config.calc_vignette_list(focus, pattern, theta, r)
        scale = ag_config.calc_scale(focus)

        def mm2pix(tup):
            angle, mm = tup
            radius_px = image.calc_radius_xy(ctr_x, ctr_y, mm * scale)
            if focus in ('CS', 'CS_IR', 'CS_OPT'):
                angle += pa_deg
            elif focus in ('NS_IR', 'NS_OPT',):
                angle += pa_deg * 2.0

            rad = math.pi * angle / 180.0
            if focus in ('CS', 'CS_IR', 'CS_OPT', 'NS_OPT'):
                x1 = int(ctr_x + radius_px * math.cos(rad))
            elif focus in ('NS_IR', ):
                x1 = int(ctr_x - radius_px * math.cos(rad))
            else:
                raise VGWError("Focus '%s' not yet implemented!" % (
                    focus))

            y1 = int(ctr_y - radius_px * math.sin(rad))
            return (x1, y1)

        return list(map(mm2pix, vlist))

    def draw(self, pluginObj):
        canvas = pluginObj.get_canvas()
        p = self.p
        thickness = self.ring_thickness
        color = self.colors

        # Draw AG probe outer movable area fov
        canvas.add(
            cvtypes.Circle(p.ctr_x, p.ctr_y, self.outer_radius,
                           color=color.outer,
                           linestyle='dash', linewidth=thickness),
            redraw=False)

        ## # Draw AG probe inner movable area fov
        ## canvas.add(
        ##     cvtypes.Circle(p.ctr_x, p.ctr_y, self.inner_radius,
        ##                        color=color.inner,
        ##                        linestyle='dash', linewidth=thickness),
        ##     redraw=False)

        # Draw AG probe position as a circle
        canvas.add(
            cvtypes.Circle(p.probe_x, p.probe_y, self.probe_radius,
                           color=color.probe,
                           linestyle='dash', linewidth=thickness),
            redraw=False)

        # Draw vignette map
        self.vig_obj = cvtypes.Polygon(self.vignette_map,
                                       color=color.vignette,
                                       linestyle='dash',
                                       linewidth=thickness)
        canvas.add(self.vig_obj, redraw=False)


class GENERICfov(TELESCOPEfov):
    """A Generic instrument FOV object.
    """

    def __init__(self, logger, image, p):
        super(GENERICfov, self).__init__(logger, image, p)

        # calculate radius of instrument fov
        inst_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, p.inst_fov)
        self.logger.debug("Instrument radius is %d pixels." % (
            inst_radius))

        self.inst_radius = inst_radius


    def draw(self, pluginObj):
        super(GENERICfov, self).draw(pluginObj)
        canvas = pluginObj.get_canvas()

        p = self.p
        thickness = self.ring_thickness
        color = self.colors

        # Draw instrument fov
        canvas.add(
            cvtypes.Circle(p.ctr_x, p.ctr_y, self.inst_radius,
                           color=color.inst,
                           linestyle='dash', linewidth=thickness),
            redraw=False)


class SHfov(object):
    """A Generic Shack-Hartmann FOV object.
    """

    def __init__(self, logger, image, p):
        self.logger = logger
        self.p = p

        # thickness of the marked rings for FOV
        self.ring_thickness = 2
        self.colors = Bunch.Bunch(inst='magenta', outer='red', inner='red',
                                  vignette='green', probe='cyan')


    def draw(self, pluginObj):
        canvas = pluginObj.get_canvas()
        p = self.p
        thickness = self.ring_thickness
        color = self.colors

        # Draw prove movable area fov
        canvas.add(
            cvtypes.Circle(p.ctr_x, p.ctr_y, p.outer_radius,
                               color=color.outer,
                               linestyle='dash', linewidth=thickness),
            redraw=False)

        # Draw coincentric catalog radii every 0.5 deg
        for cat_radius in p.cat_radii:
            canvas.add(
                cvtypes.Circle(p.ctr_x, p.ctr_y, cat_radius,
                                   color='white',
                                   linestyle='dash', linewidth=thickness),
                redraw=False)


class MOIRCSfov(TELESCOPEfov):
    """A MOIRCS (instrument) FOV object.
    """

    def __init__(self, logger, image, p):
        super(MOIRCSfov, self).__init__(logger, image, p)

        # TODO: read these from a config file
        fov = 0.033000
        # moircs fov 7arcmin x 4 arcmin(0.1166 x 0.0666 in degree) half values
        fov_hh = 0.0583333
        fov_hw = 0.0333333
        # moircs fov with vignette  9.2 arcmin x 6.2 arcmin  half values
        vig_hh = 0.0766666
        vig_hw = 0.05166667

        # Figure out rotation of MOIRCS rectangular FOV and chip markings
        # NOTE: ginga plots by WCS, so we don't need to account for
        # rotation of the DSS image for now
        ## header = image.get_header()
        ## #dssrot_deg = image.get_wcs_rotation_deg()
        ## ((xrot_ref, yrot_ref),
        ##  (cdelt1_ref, cdelt2_ref)) = wcs.get_xy_rotation_and_scale(header)
        ## dssrot_deg = yrot_ref
        ## self.logger.info("Image rotation=%f, pa=%f cdelt1=%f" % (
        ##         dssrot_deg, p.ag_pa, cdelt1_ref))

        # Angle we should draw the object at is therefore
        theta = p.ag_pa
        #self.theta = -theta
        self.theta = theta
        self.logger.debug("rotation is %f deg" % (self.theta))

        # coords of the MOIRCS FOV
        self.fov_pts = []
        for (xoff, yoff) in ((-fov_hw, -fov_hh), (-fov_hw, +fov_hh),
                             (+fov_hw, +fov_hh), (+fov_hw, -fov_hh)):
            self.fov_pts.append(image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                    xoff, yoff))
        # coords of FOV w/vignette
        self.vig_pts = []
        for (xoff, yoff) in ((-vig_hw, -vig_hh), (-vig_hw, +vig_hh),
                             (+vig_hw, +vig_hh), (+vig_hw, -vig_hh)):
            self.vig_pts.append(image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                    xoff, yoff))

        # calculate the position of the MOIRCS CCD chip and indicate
        # the position by text
        self.c1x, self.c1y = image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                 0, -fov_hh/2.0)
        self.c2x, self.c2y = image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                 0, +fov_hh/2.0)


    def draw(self, pluginObj):
        super(MOIRCSfov, self).draw(pluginObj)
        canvas = pluginObj.get_canvas()

        thickness = self.ring_thickness
        color = self.colors

        obj = cvtypes.CompoundObject()
        canvas.add(obj, redraw=False)

        obj.add_object(cvtypes.Polygon(self.fov_pts,
                                       color=color.inst,
                                       linestyle='dash',
                                       linewidth=thickness))
        obj.add_object(cvtypes.Polygon(self.vig_pts,
                                       color='blue',
                                       linestyle='solid',
                                       linewidth=thickness))
        obj.add_object(cvtypes.Text(self.c1x, self.c1y, "Chip1",
                                    color='white'))
        obj.add_object(cvtypes.Text(self.c2x, self.c2y, "Chip2",
                                    color='white'))
        obj.rotate(self.theta, xoff=self.p.ctr_x, yoff=self.p.ctr_y)

        # MOIRCS is offset by 45 deg wrt cassegrain flange.  Standard
        # vignette map needs to be rotated to be properly aligned
        vig_rot = 2.0 * self.theta - 45.0
        self.vig_obj.rotate(vig_rot, xoff=self.p.ctr_x, yoff=self.p.ctr_y)


class SPCAMfov(object):
    """A SPCAM (instrument) FOV object.
    """

    def __init__(self, logger, image, p):
        self.logger = logger
        self.p = p
        outer_fov = 0.2879
        ## inner_fov = 0.004167
        #P_OPT   FOV 0.46    -0.018702048    0.47    -0.072          0.072
        fov = 0.46
        pb_minus = 0.018702048
        pb_plus  = 0.47
        pb_height = 0.072

        # thickness of the marked rings for FOV
        self.ring_thickness = 2
        self.colors = Bunch.Bunch(inst='magenta', outer='red', inner='red',
                                  vignette='green', probe='cyan')

        ## # calculate radius of probe inner movable area fov
        ## inner_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, inner_fov)
        ## self.logger.debug("Probe inner movable area radius is %d pixels." % (
        ##     inner_radius))

        # calculate probe circle
        probe_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, p.probe_head_fov)
        self.logger.debug("Actual probe head radius is %d pixels." % (
            probe_radius))

        # calculate radius of probe outer movable area fov
        outer_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, outer_fov)
        self.logger.debug("Probe outer movable area radius is %d pixels." % (
            outer_radius))

        # calculate radius of instrument fov
        inst_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, p.inst_fov)
        self.logger.debug("Instrument radius is %d pixels." % (
            inst_radius))

        self.inst_radius = inst_radius
        ## self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.probe_radius = probe_radius

        # Figure out rotation of SPCAM drawings
        # NOTE: ginga plots by WCS, so we don't need to account for
        # rotation of the DSS image for now
        ## header = image.get_header()
        ## ((xrot_ref, yrot_ref),
        ##  (cdelt1_ref, cdelt2_ref)) = wcs.get_xy_rotation_and_scale(header)
        ## dssrot_deg = yrot_ref

        theta = p.ag_pa
        # change direction of rotation for P_OPT
        self.theta = -theta
        self.logger.debug("rotation is %f deg" % (self.theta))

        # coords of the MOVEABLE PROBE FOV
        self.fov_pts = []
        for (xoff, yoff) in ((-pb_minus, -pb_height), (-pb_minus, +pb_height),
                             (+pb_plus, +pb_height), (+pb_plus, -pb_height)):
            self.fov_pts.append(image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                    xoff, yoff))

    def draw(self, pluginObj):
        self.pluginObj = pluginObj
        canvas = pluginObj.get_canvas()
        thickness = self.ring_thickness
        color = self.colors
        p = self.p

        obj = cvtypes.CompoundObject()
        canvas.add(obj, redraw=False)

        # Draw instrument fov
        obj.add_object(
            cvtypes.Circle(p.ctr_x, p.ctr_y, self.inst_radius,
                           color=color.inst,
                           linestyle='dash', linewidth=thickness))

        # Draw outer telescope foci fov
        obj.add_object(
            cvtypes.Circle(p.ctr_x, p.ctr_y, self.outer_radius,
                           color=color.outer,
                           linestyle='dash', linewidth=thickness))

        # Draw probe position
        obj.add_object(
            cvtypes.Circle(p.probe_x, p.probe_y, self.probe_radius,
                           color=color.probe,
                           linestyle='dash', linewidth=thickness))
        ## # Draw AG probe inner movable area fov
        ## obj.addObject(
        ##     cvtypes.Circle(p.ctr_x, p.ctr_y, self.inner_radius,
        ##                        color=color.inner,
        ##                        linestyle='dash', linewidth=thickness))

        # Draw probe movable area
        self.prb_area = cvtypes.Polygon(self.fov_pts,
                                        color=color.probe,
                                        linestyle='dash',
                                        linewidth=thickness)
        obj.add_object(self.prb_area)

        # Rotate according to PA
        obj.rotate(self.theta, xoff=self.p.ctr_x, yoff=self.p.ctr_y)

    def filter_results(self, starlist):
        return self.pluginObj.filter_results(starlist, self.prb_area)


class HSCfov(object):
    """A Hyper Suprime-Cam (instrument) FOV object.
    """

    def __init__(self, logger, image, p):
        self.logger = logger
        self.p = p
        self.ccds_obj = None
        self.dith_obj = None

        # thickness of the marked rings for FOV
        self.ring_thickness = 2
        self.colors = Bunch.Bunch(inst='magenta', outer='red', inner='red',
                                  vignette='green', probe='cyan')


    def draw(self, pluginObj):
        canvas = pluginObj.get_canvas()
        self.pluginObj = pluginObj
        p = self.p
        thickness = self.ring_thickness
        color = self.colors

        gons = []
        for points in p.polygons:
            # Add CCD image polygon
            self.logger.info("Plotting polygon %s" % str(points))
            gons.append(cvtypes.Polygon(points,
                                        color=color.inst,
                                        linestyle='dash',
                                        linewidth=thickness))
        self.ccds_obj = cvtypes.CompoundObject(*gons)
        canvas.add(self.ccds_obj, redraw=False)

        pixgons = []
        for points in p.agarea_pixel_polygons:
            # Add internal dither image polygon
            self.logger.info("Plotting dithering polygon %s" % str(points))
            pixgons.append(cvtypes.Polygon(points,
                                           color=color.vignette,
                                           linestyle='dash',
                                           linewidth=1))
        self.dith_obj = cvtypes.CompoundObject(*pixgons)
        canvas.add(self.dith_obj, redraw=False)

        # for (x, y, r) in p.circles:
        #     canvas.add(
        #         cvtypes.Circle(x, y, r,
        #                            color=color.probe,
        #                            linestyle='solid', linewidth=1),
        #         redraw=False)

    def filter_results(self, starlist):
        starlist = self.pluginObj.filter_results(starlist, self.dith_obj)

        self.logger.debug('all STARLIST=%s' % (str(starlist)))

        p = self.p
        return self.hsc_filter_candidates(p, p.queries, starlist)

    def hsc_filter_candidates(self, p, queries, all_stars):
        # interesting items:
        # p.ag_pa, p.limitmag, p.goodmag
        # p.polygons, p.circles
        # queries: ((ra, dec, radius), ... ) 1 for each ccd
        # all_stars: query result for each circle

        fabs = math.fabs

        goodmag = p['goodmag']
        limitmag = p['limitmag']
        bright_end = 1.0
        too_bright = 10.0
        best_flag = 2.0

        # TEMP: until we fix the query
        #all_stars = filter(lambda star: star['mag'] <= limitmag, all_stars)

        for star in all_stars:

            pref = 0

            diff_mag = goodmag - star['mag']

            if diff_mag > bright_end:
                pref +=  too_bright
            else:
                pref += fabs(diff_mag)

            pref += fabs(best_flag-star['flag'])

            star['preference'] = pref
            if (star['description'] is not None and
                'BLACKLISTED' in star['description']):
                star['preference'] = 9999999

        all_stars = sorted(all_stars, key=itemgetter('preference'))
        reset_priority(all_stars)

        return all_stars


class MIMIZUKUfov(TELESCOPEfov):
    """A MIMIZUKU (instrument) FOV object.

    Kamizuka: We confirmed that the field stacker (FS) rotator rotates in
    the same direction as the instrument rotator. Then, they are additive.

    But, their rotation direction is opposite to the definition of the
    position angle. Then, the position angle of the prohibited area will
    be -[InR angle + FS angle].

    In addition to those, an offset will be added, but this offset has an
    uncertainty.  It depends on the InR angle in the instrument exchange.
    Okita-san says that it may be -90 deg.
    """

    def __init__(self, logger, image, p):
        super(MIMIZUKUfov, self).__init__(logger, image, p)

        # offset in attachment to the Cassegrain flange
        self.mmz_offset = -90.0
        # TODO: read these from a config file
        fov = 0.033000
        # mimizuku fov with detector vignette  9.2 arcmin x 2.4 arcmin
        #  (half values)
        vig_hh = 0.02
        vig_hw = 0.0766666

        # Angle we should draw the object at is therefore
        # TODO: need to account for MIMIZUKU field stacker rotation
        #   plus offset
        theta = p.ag_pa
        self.theta = -theta
        self.logger.debug("rotation is %f deg" % (self.theta))

        # coords of detector vignette
        self.vig_pts = []
        for (xoff, yoff) in ((-vig_hw, -vig_hh), (-vig_hw, +vig_hh),
                             (+vig_hw, +vig_hh), (+vig_hw, -vig_hh)):
            self.vig_pts.append(image.add_offset_xy(p.ctr_x, p.ctr_y,
                                                    xoff, yoff))

        # calculate the position of the MIMIZUKU CCD chip and indicate
        # the position by text
        ## self.c1x, self.c1y = image.add_offset_xy(p.ctr_x, p.ctr_y,
        ##                                          0, -fov_hh/2.0)
        ## self.c2x, self.c2y = image.add_offset_xy(p.ctr_x, p.ctr_y,
        ##                                          0, +fov_hh/2.0)

        self.dith_obj = None

    def draw(self, pluginObj):
        super(MIMIZUKUfov, self).draw(pluginObj)
        canvas = pluginObj.get_canvas()
        self.pluginObj = pluginObj

        thickness = self.ring_thickness
        color = self.colors

        obj = cvtypes.CompoundObject()
        canvas.add(obj, redraw=False)

        self.dith_obj = cvtypes.Polygon(self.vig_pts,
                                        color='blue',
                                        linestyle='solid',
                                        linewidth=thickness)
        obj.add_object(self.dith_obj)
        # obj.add_object(cvtypes.Text(self.c1x, self.c1y, "Chip1",
        #                             color='white'))
        # obj.add_object(cvtypes.Text(self.c2x, self.c2y, "Chip2",
        #                             color='white'))
        obj.rotate(self.theta, xoff=self.p.ctr_x, yoff=self.p.ctr_y)

        # MIMIZUKU is offset wrt cassegrain flange.  Standard
        # vignette map needs to be rotated to be properly aligned
        vig_rot = 2.0 * self.theta - self.mmz_offset
        self.vig_obj.rotate(vig_rot, xoff=self.p.ctr_x, yoff=self.p.ctr_y)

    def filter_results(self, starlist):
        prohibited = self.pluginObj.filter_results(starlist, self.dith_obj)

        good = set(starlist) - set(prohibited)
        self.logger.info("len prohibited=%d  good=%d" % (
            len(prohibited), len(good)))

        starlist = list(sorted(good, key=itemgetter('preference'),
                                  reverse=True))

        # prohibited stars get added after
        prohibited = list(sorted(prohibited, key=itemgetter('preference'),
                                  reverse=True))
        starlist.extend(prohibited)

        # finally, renumber priorities as position in list
        reset_priority(starlist)

        self.logger.debug('all STARLIST=%s' % (str(starlist)))
        return starlist


def reset_priority(updated_stars):
    for num, star in enumerate(updated_stars):
        num += 1
        if 'priority' not in star or (star['priority'] is None):
            star['priority'] = num


#END
