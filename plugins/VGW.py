#
# VGW.py -- VGW plugin for fits viewer
# 
# Eric Jeschke (eric@naoj.org)
#
import math
import time
import gtk
import os
import numpy
from operator import itemgetter

# Ginga imports
from ginga import GingaPlugin
from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
from ginga.misc import Future, Bunch
from ginga import AstroImage
from ginga.util import wcs, dp


# PYHOME imports
import remoteObjects as ro
import astro.radec as radec
# For HSC AG drawing
import astro.wcs as astro_wcs
import cfg.g2soss as g2soss
import SOSS.GuiderInt.ag_config as ag_config

# Local application imports
from util import AGCCDPositions
from util import g2calc

class VGWError(Exception):
    pass

msg_auto_success = """Automatic guide star selection succeeded."""
msg_semiauto_failure = """Automatic guide star selection failed!

Please select a guide star manually."""
msg_auto_manual = """Manual mode selected:

Please select a guide star manually."""

# Where sounds are stored
soundhome = os.path.join(g2soss.soundhome, 'ogg', 'en')

# Sounds
snd_auto_failure = os.path.join(soundhome, "auto_guide_star_sel_failed5.ogg")
snd_auto_select_manual = os.path.join(soundhome,
                                      "select_guide_star_manually7.ogg")

snd_region_failure = os.path.join(soundhome, "auto_regionselect_failed4.ogg")
snd_region_select_manual = os.path.join(soundhome,
                                      "select_region_manually3.ogg")

snd_agarea_failure = os.path.join(soundhome, "auto_agareaselect_failed3.ogg")
snd_agarea_select_manual = os.path.join(soundhome,
                                      "select_area_manually3.ogg")


class VGW(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(VGW, self).__init__(fv)

        self.count = 0

        self.colorcalc = 'cyan'
        self.qdaschname = 'QDAS_VGW'

        self.catalogs = self.fv.get_ServerBank()

        # load preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.createCategory('plugin_VGW')
        self.settings.load(onError='silent')

        self.dss_server_sfx = self.settings.get('dss_server_suffix',
                                                '@subaru')
        
        # For drawing HSC guide CCDs
        self.hscguideinfo = AGCCDPositions.SCAGCCDPositions()

        # For image FWHM type calculations
        self.iqcalc = g2calc.IQCalc(self.logger)
        
        # For the SV_DRIVE command
        self.sv_dst_x = None
        self.sv_dst_y = None
        self.sv_obj_x = None
        self.sv_obj_y = None

    #############################################################
    #    Here come the VGW commands
    #############################################################


    def mark_position(self, tag, future,
                      v_lan_data=None, x=None, y=None, mode=None,
                      mark=None, size=None, color=None):
        # X and Y are in FITS data coordinates
        x, y = x-1, y-1

        chname = v_lan_data.upper()
        chinfo = self.fv.get_channelInfo(chname)

        p = future.get_data()

        canvas = chinfo.fitsimage
        try:
            self._mark(chname, canvas, x, y, mode, mark, size, color)
            p.setvals(result='ok')

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))
            
        future.resolve(0)
        
                
    def region_selection(self, tag, future, motor=None, v_lan_data=None,
                         select_mode=None, x_region=None, y_region=None,
                         exptime=3000):
    

        # Copy the current ag image specified by (v_lan_data) into the
        # destination channel
        chname = v_lan_data
        if not self.fv.has_channel(self.qdaschname):
            self.fv.add_channel(self.qdaschname)

        chinfo = self.fv.get_channelInfo(self.qdaschname)

        # Deactivate plugin if one is already running
        pluginName = 'Region_Selection'
        if chinfo.opmon.is_active(pluginName):
            chinfo.opmon.deactivate(pluginName)
            self.fv.update_pending()

        self.copy_data(chname, self.qdaschname)
        self.fv.ds.raise_tab(self.qdaschname)

        image = chinfo.fitsimage.get_image()
        assert image is not None, \
               VGWError("Null image for '%s'!" % (self.qdaschname))

        try:
            exptime = int(exptime)
        except:
            exptime = 3000

        p = future.get_data()
        p.setvals(exptime=exptime, auto=False, recenter=True,
                  obj_x=0.0, obj_y=0.0, fwhm=0.0,
                  skylevel=0.0, brightness=0.0,
                  dx=x_region // 2, dy=y_region // 2)

        # TODO: get old values
        #x, y, exptime, fwhm, brightness, skylevel, objx, objy
        ## p = self.get_default_parameters('region_selection',
        ##                             ('x', 'y', 'exptime', 'fwhm', 'brightness',
        ##                             'skylevel', 'objx', 'objy'))

        p.x, p.y = image.width // 2, image.height // 2
        p.obj_x, p.obj_y = p.x, p.y
        p.x1, p.y1 = max(0, p.x - p.dx), max(0, p.y - p.dy)
        p.x2 = min(image.width-1,  p.x + p.dx)
        p.y2 = min(image.height-1, p.y + p.dy)

        rsinfo = chinfo.opmon.getPluginInfo(pluginName)
        rsobj = rsinfo.obj

        thr = rsobj.threshold
        if (thr == ''):
            thr = None
        p.setvals(alg='v1', radius=rsobj.radius, threshold=thr)

        select_mode = select_mode.lower()
        if select_mode in ('auto', 'semiauto'):
            try:
                self._auto_region_selection(image, p)

                rsobj.show_region(p)
                p.exptime = rsobj.exptime

                self.map_back_to_ccd(p, image)
                future.resolve(0)
                return

            except Exception, e:
                errmsg = "Automatic region selection failed: %s" % str(e)
                self.logger.error(errmsg)
                self.fv.play_soundfile(snd_region_failure, priority=19)

        #elif select_mode in ('manual', ):

        # Need to create a second future here only because we need to
        # do a map_back_to_ccd with the result before returning results...
        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._region_selection_cb, future,
                             p, image)
        
        # Invoke the operation manually
        chinfo.opmon.start_plugin_future(self.qdaschname, pluginName,
                                         future2)
        self.fv.update_pending(timeout=0.10)
        self.fv.play_soundfile(snd_region_select_manual, priority=20)


    def _region_selection_cb(self, future2, future, p, image):
        try:
            if p.result == 'ok':
                self.map_back_to_ccd(p, image)

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        self.logger.debug("region selection cb terminating: res=%s" % (str(p)))
        future.resolve(0)
    
    def _auto_region_selection(self, image, p):

        x1, y1 = 0, 0
        x2, y2 = image.width-1, image.height-1

        if p.alg == 'v2':
            qualsize = self.iqcalc.qualsize
        else:
            # NOTE: qualsize_old is Kosugi-san's old algorithm
            qualsize = self.iqcalc.qualsize_old

        qs = qualsize(image, x1, y1, x2, y2,
                      radius=p.radius, threshold=p.threshold)
        p.x, p.y = qs.x, qs.y

        # set region from center of detected object
        p.x1, p.y1 = max(0, p.x - p.dx), max(0, p.y - p.dy)
        p.x2 = min(image.width-1,  p.x + p.dx)
        p.y2 = min(image.height-1, p.y + p.dy)

        p.fwhm = qs.fwhm
        p.brightness = qs.brightness
        p.skylevel = qs.skylevel
        p.obj_x = qs.objx
        p.obj_y = qs.objy
        p.result = 'ok'

        self.logger.info("Region selection succeeded %s" % (
            str(p)))


    def ag_guide_area_selection(self, tag, future,
                                motor=None, v_lan_data=None,
                                select_mode=None, x_region=None, y_region=None,
                                ag_area=None):
    
        v_lan_data = v_lan_data.upper()

        # Copy the current ag image specified by (v_lan_data) into the
        # destination channel
        chname = v_lan_data
        if not self.fv.has_channel(self.qdaschname):
            self.fv.add_channel(self.qdaschname)

        chinfo = self.fv.get_channelInfo(self.qdaschname)

        # Deactivate plugin if one is already running
        pluginName = 'AgAreaSelection'
        if chinfo.opmon.is_active(pluginName):
            chinfo.opmon.deactivate(pluginName)
            self.fv.update_pending()

        self.copy_data(chname, self.qdaschname)
        self.fv.ds.raise_tab(self.qdaschname)

        image = chinfo.fitsimage.get_image()
        assert image is not None, \
               VGWError("Null image for '%s'!" % (self.qdaschname))

        ag_area = ag_area.lower()
        
        # Get the exposure area from the AG header
        agh = Bunch.Bunch(image.get('agheader'))
        self.logger.info("AG header is %s" % (str(agh)))
        if agh.has_key('expRangeX'):
            # Need to convert via ccd2pix here?
            er_x1 = agh.expRangeX
            er_y1 = agh.expRangeY
            er_x2 = er_x1 + agh.expRangeDX - 1
            er_y2 = er_y1 + agh.expRangeDY - 1
        
        p = future.get_data()
        p.setvals(exptime=0.0, auto=False,
                  obj_x=0.0, obj_y=0.0, fwhm=0.0,
                  skylevel=0.0, brightness=0.0,
                  ag_area=ag_area, agkey=v_lan_data,
                  er_x1=er_x1, er_y1=er_y1,
                  er_x2=er_x2, er_y2=er_y2,
                  dx=x_region // 2, dy = y_region // 2)

        # TODO: get old values
        #x, y, exptime, fwhm, brightness, skylevel, objx, objy
        ## p = self.get_default_parameters('region_selection',
        ##                             ('x', 'y', 'exptime', 'fwhm', 'brightness',
        ##                             'skylevel', 'objx', 'objy'))

        p.x, p.y = image.width // 2, image.height // 2
        p.x1, p.y1 = max(0, p.x - p.dx), max(0, p.y - p.dy)
        p.x2 = min(image.width-1,  p.x + p.dx)
        p.y2 = min(image.height-1, p.y + p.dy)

        rsinfo = chinfo.opmon.getPluginInfo(pluginName)
        rsobj = rsinfo.obj

        thr = rsobj.threshold
        if (thr == ''):
            thr = None
        p.setvals(alg='v1', radius=rsobj.radius, threshold=thr)

        select_mode = select_mode.lower()
        if select_mode in ('auto', 'semiauto'):
            try:
                self._auto_region_selection(image, p)

                rsobj.show_area(p)
                p.exptime = rsobj.exptime

                self.map_back_to_ccd(p, image)
                future.resolve(0)
                return

            except Exception, e:
                errmsg = "Automatic ag area selection failed: %s" % str(e)
                self.logger.error(errmsg)
                self.fv.play_soundfile(snd_agarea_failure, priority=19)

        #elif select_mode in ('manual', ):
            
        # Need to create a second future here only because we need to
        # do a map_back_to_ccd with the result before returning results...
        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._agarea_selection_cb, future,
                             p, image)
        
        # Invoke the operation manually
        chinfo.opmon.start_plugin_future(self.qdaschname, pluginName, future2)
        self.fv.update_pending(timeout=0.10)
        self.fv.play_soundfile(snd_agarea_select_manual, priority=20)


    def _agarea_selection_cb(self, future2, future, p, image):
        try:
            if p.result == 'ok':
                self.map_back_to_ccd(p, image)

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        self.logger.debug("region selection cb terminating: res=%s" % (str(p)))
        future.resolve(0)
        

    def sv_drive(self, tag, future,
                 motor=None, slit_x=None, slit_y=None,
                 object_x=None, object_y=None):

        chname = 'SV'
        if not self.fv.has_channel(self.qdaschname):
            self.fv.add_channel(self.qdaschname)

        # Invoke the operation manually
        chinfo = self.fv.get_channelInfo(self.qdaschname)

        # Deactivate plugin if one is already running
        pluginName = 'Sv_Drive'
        if chinfo.opmon.is_active(pluginName):
            chinfo.opmon.deactivate(pluginName)
            self.fv.update_pending()

        self.copy_data(chname, self.qdaschname)
        self.fv.ds.raise_tab(self.qdaschname)
        
        image = chinfo.fitsimage.get_image()
        assert image is not None, \
               VGWError("Null image for '%s'!" % (self.qdaschname))

        agh = Bunch.Bunch(image.get('agheader'))
        
        # Convert pixel coords on image back to CCD coords
        obj_x, obj_y = None, None
        dst_x, dst_y = None, None
        if self.sv_obj_x is not None:
            obj_x, obj_y = self.ccd2pix(self.sv_obj_x, self.sv_obj_y, agh.binX,
                                        agh.expRangeX, agh.expRangeY)
            dst_x, dst_y = self.ccd2pix(self.sv_dst_x, self.sv_dst_y, agh.binX,
                                        agh.expRangeX, agh.expRangeY)
        
        # Set defaults and adjust for difference between data coords and
        # fits coords
        if slit_x is not None:
            slit_x -= 1
        else:
            slit_x = dst_x
        if slit_y is not None:
            slit_y -= 1
        else:
            slit_y = dst_y
        if object_x is not None:
            object_x -= 1
        else:
            object_x = obj_x
        if object_y is not None:
            object_y -= 1
        else:
            object_y = obj_y

        p = future.get_data()
        p.setvals(dst_x=slit_x, dst_y=slit_y,
                  obj_x=object_x, obj_y=object_y,
                  x1=0, y1=0, x2=image.width-1, y2=image.height-1)

        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._sv_drive_cb, future,
                             p, image)
        
        chinfo.opmon.start_plugin_future(self.qdaschname, pluginName,
                                         future2)

    def _sv_drive_cb(self, future2, future, p, image):

        try:
            if p.result == 'ok':
                # save positions mapped back to AG CCD for future
                # incantations of the same plugin
                agh = Bunch.Bunch(image.get('agheader'))
                self.sv_obj_x, self.sv_obj_y = self.pix2ccd(p.obj_x, p.obj_y,
                                 agh.binX, agh.expRangeX, agh.expRangeY)
                self.sv_dst_x, self.sv_dst_y = self.pix2ccd(p.dst_x, p.dst_y,
                                 agh.binX, agh.expRangeX, agh.expRangeY)

                dst_ra_deg, dst_dec_deg = image.pixtoradec(p.dst_x, p.dst_y)
                obj_ra_deg, obj_dec_deg = image.pixtoradec(p.obj_x, p.obj_y)

                p.dst_x += 1
                p.dst_y += 1
                p.equinox = 2000.0
                p.dst_equinox = 2000.0
                p.obj_equinox = 2000.0
                p.dst_ra = radec.raDegToString(dst_ra_deg)
                p.dst_dec = radec.decDegToString(dst_dec_deg)
                p.obj_x += 1
                p.obj_y += 1
                p.obj_ra = radec.raDegToString(obj_ra_deg)
                p.obj_dec = radec.decDegToString(obj_dec_deg)

                sep_ra, sep_dec = wcs.get_RaDecOffsets(obj_ra_deg, obj_dec_deg,
                                                       dst_ra_deg, dst_dec_deg)
                p.rel_ra = radec.offsetRaDegToString(sep_ra)
                p.rel_dec = radec.decDegToString(sep_dec)
        
        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        self.logger.debug("sv_drive cb terminating: res=%s" % (str(p)))
        future.resolve(0)
    

    def ag_auto_select(self, tag, future,
                       motor=None, equinox=None, ra=None, dec=None,
                       f_select=None, probe_ra=None, probe_dec=None,
                       ag_pa=None, probe_r=None, probe_theta=None,
                       probe_x=None, probe_y=None, select_mode=None,
                       dss_mode=None, ag_area=None, instrument_name=None,
                       limitmag=None, goodmag=None, fov_pattern=None):

        chname = 'DSS'
        if not self.fv.has_channel(chname):
            self.fv.add_channel(chname)
        chinfo = self.fv.get_channelInfo(chname)

        ra_deg = radec.funkyHMStoDeg(ra)
        dec_deg = radec.funkyDMStoDeg(dec)

        # Calculate fov from instrument
        magic_constant = 2.5
        outer_fov = ag_config.calc_outer_fov(f_select)
        inst_fov = ag_config.calc_instrument_fov(instrument_name)
        probe_head_fov = ag_config.calc_probe_head_fov(f_select)
        probe_vignette_fov = ag_config.calc_probe_vignette_fov(f_select)
        
        # For DSS, ra and dec are specified in traditional format
        ra_txt = radec.raDegToString(ra_deg, format='%02d:%02d:%06.3f')
        dec_txt = radec.decDegToString(dec_deg,
                                       format='%s%02d:%02d:%05.2f')

        if not dss_mode:
            dss_mode = 'off'
        else:
            dss_mode = dss_mode.lower()
            if dss_mode == 'on':
                dss_mode = 'dss1'
                
        # DSS FOV (deg)
        dss_fov_deg = outer_fov * magic_constant
        # Catalog FOV (deg)
        cat_fov_deg = outer_fov

        p = future.get_data()
        p.setvals(image=None, ra_deg=ra_deg, dec_deg=dec_deg,
                  inst_fov=inst_fov, outer_fov=outer_fov,
                  probe_head_fov=probe_head_fov,
                  probe_vignette_fov=probe_vignette_fov,
                  cat_fov_deg=cat_fov_deg,
                  ag_pa=ag_pa, f_select=f_select,
                  select_mode=select_mode, probe_ra=probe_ra,
                  probe_dec=probe_dec, probe_r=probe_r, probe_theta=probe_theta,
                  probe_x=probe_x, probe_y=probe_y, dss_mode=dss_mode,
                  ag_area=ag_area, instrument_name=instrument_name,
                  fov_pattern=fov_pattern, limitmag=limitmag, goodmag=goodmag,
                  chname=chname, equinox=equinox)
        
        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._ag_auto_select_cont3, future)
        
        # Open up the UI
        if not chinfo.opmon.is_active('AgAutoSelect'):
            chinfo.opmon.start_plugin_future(chname, 'AgAutoSelect',
                                             future2, alreadyOpenOk=True)
        else:
            self.fv.ds.raise_tab('DSS')

        def get_dss_image(p):
            pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
            pluginObj = pluginInfo.obj

            # Assume square image?
            wd_deg = dss_fov_deg
            ht_deg = dss_fov_deg

            # width and height are specified in arcmin
            sgn, deg, mn, sec = radec.degToDms(wd_deg)
            wd = deg*60.0 + float(mn) + sec/60.0
            sgn, deg, mn, sec = radec.degToDms(ht_deg)
            ht = deg*60.0 + float(mn) + sec/60.0
            ## sgn, deg, mn, sec = radec.degToDms(radius_deg)
            ## radius = deg*60.0 + float(mn) + sec/60.0

            # Get the DSS image into the window
            try:
                if dss_mode != 'off':
                    params = dict(ra=ra_txt, dec=dec_txt, width=wd, height=ht)

                    # Get the DSS server name
                    default_server = dss_mode + self.dss_server_sfx
                    server = self.settings.get('dss_server', default_server)

                    # Query the server and download file
                    fitspath = pluginObj.get_sky_image(server, params)
                    image = self.fv.load_file(fitspath, chname=chinfo.name)
                else:
                    # make blank image
                    #px_scale = self.calculate_fov_scale(p.f_select)
                    px_scale = 0.0004722298
                    #px_scale = 0.000280178318866
                    image = dp.create_blank_image(ra_deg, dec_deg,
                                                  dss_fov_deg,
                                                  px_scale, 0.0,
                                                  cdbase=[-1, 1],
                                                  logger=self.logger)
                    self.fv.gui_call(chinfo.fitsimage.set_image, image)

                p.image = image

            except Exception as e:
                errmsg = "Error querying dss server: %s" % (str(e))
                p.image = None
                self.logger.error(errmsg)
                self.fv.show_error(errmsg)

        # Offload this network task to a non-gui thread
        f_dss = Future.Future()
        f_dss.freeze(get_dss_image, p)

        # Clear old data from canvas
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        pluginObj.reset()
        
        chinfo.fitsimage.onscreen_message("Querying image db...",
                                          delay=1.0)
        self.fv.showStatus("Querying sky image database ...")
        f_rest = Future.Future()
        f_rest.freeze(self.fv.gui_do, self._ag_auto_select_cont1, future2,
                      future)
        f_dss.add_callback('resolved', f_rest.thaw)
        self.fv.nongui_do_future(f_dss)

    def _ag_auto_select_cont1(self, future2, future):
        self.logger.debug("continuation 1 resumed...")
        p = future.get_data()
        if p.image is None:
            # TODO: pop up an error message
            self.fv.show_error("No DSS image returned!")
            future.resolve(-1)
            return
        self.fv.showStatus("Got DSS image.")

        image = p.image
        # Now that we have an image, we can do some WCS calculations
        
        probe_ra_deg = radec.funkyHMStoDeg(p.probe_ra)
        probe_dec_deg = radec.funkyDMStoDeg(p.probe_dec)
        probe_x, probe_y = image.radectopix(probe_ra_deg, probe_dec_deg)

        # calculate center pixel for ra/dec
        ctr_x, ctr_y = image.radectopix(p.ra_deg, p.dec_deg)
        # calculate radius of probe vignetting
        probe_vignette_radius = image.calc_radius_xy(ctr_x, ctr_y,
                                                     p.probe_vignette_fov)

        p.setvals(ctr_x=ctr_x, ctr_y=ctr_y, probe_x=probe_x, probe_y=probe_y,
                  probe_vignette_radius=probe_vignette_radius,
                  probe_ra_deg=probe_ra_deg, probe_dec_deg=probe_dec_deg)

        # Catalog function doesn't like to see CS_OPT or CS_IR, just "CS"
        f_select0 = p.f_select.upper()
        if f_select0 in ('CS_OPT', 'CS_IR'):
            f_select0 = 'CS'
        radius = p.cat_fov_deg * 60.0
            
        def query_catalogs(p):
            try:
                # Get preferred guide star catalog for AG
                catname = self.settings.get('AG_catalog', 'ag@subaru')
                starcat = self.catalogs.getCatalogServer(catname)

                # Query catalog
                ## query_result = starcat.search_ag(
                ##     ra_deg=p.ra_deg, dec_deg=p.dec_deg, fov_deg=p.cat_fov_deg,
                ##     probe_ra_deg=probe_ra_deg, probe_dec_deg=probe_dec_deg,
                ##     focus=f_select0, inst_name=p.instrument_name,
                ##     probe_r=p.probe_r, probe_theta=p.probe_theta,
                ##     probe_x=p.probe_x, probe_y=p.probe_y, pos_ang_deg=p.ag_pa,
                ##     upper_mag=p.limitmag, pref_mag=p.goodmag,
                ##     fov_pattern=p.fov_pattern, equinox=p.equinox,
                ##     )
                query_result = starcat.search(
                    ra=str(p.ra_deg), dec=str(p.dec_deg), equinox=p.equinox,
                    r1=str(0.0), r2=str(radius),
                    m2=str(p.limitmag), m1=str(p.goodmag),
                    focus=f_select0, pa=p.ag_pa,
                    inst_name=p.instrument_name,
                    probe_ra_deg=probe_ra_deg,
                    probe_dec_deg=probe_dec_deg,
                    probe_r=p.probe_r,
                    probe_theta=p.probe_theta,
                    probe_x=p.probe_x, probe_y=p.probe_y,
                    fov_pattern=p.fov_pattern)

                info, starlist = starcat.process_result(query_result)
                p.info = info
                self.logger.debug("info=%s" % (str(info)))
                p.starlist = starlist
                if len(starlist) > 0:
                    p.selected = [ starlist[0] ]
                else:
                    p.selected = []

            except Exception, e:
                errmsg = "Error querying star catalog: %s" % (str(e))
                self.logger.error(errmsg)
                self.fv.show_error(errmsg)
                p.setvals(info={}, starlist=None, selected=[], error=errmsg,
                          image=None)
                raise VGWError(errmsg)

        # Offload this network task to a non-gui thread
        f_cat = Future.Future()
        f_cat.freeze(query_catalogs, p)

        chinfo = self.fv.get_channelInfo(p.chname)
        chinfo.fitsimage.onscreen_message("Querying catalog db...",
                                          delay=1.0)
        self.fv.showStatus("Querying catalog for objects ...")
        f_rest = Future.Future()
        f_rest.freeze(self.fv.gui_do, self._ag_auto_select_cont2, future2,
                      future)
        f_cat.add_callback('resolved', f_rest.thaw)
        self.fv.nongui_do_future(f_cat)

    def _ag_auto_select_cont2(self, future2, future):
        self.logger.debug("continuation 2 resumed...")
        p = future.get_data()

        if p.starlist is None:
            self.fv.show_error("No catalog data returned!")
            future.resolve(-1)
            return

        self.fv.showStatus("catalog returned %d results" % (
            len(p.starlist)))

        # Plot all the data
        insname = p.instrument_name.lower()
        if insname == 'moircs':
            plotObj = MOIRCSfov(self.logger, p.image, p)
        elif insname == 'spcam':
            plotObj = SPCAMfov(self.logger, p.image, p)
        else:
            plotObj = GENERICfov(self.logger, p.image, p)
        
        chinfo = self.fv.get_channelInfo(p.chname)
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        # TEMP: fix
        #pluginObj.limit_stars_to_area = True
        pluginObj.limit_stars_to_area = False
        pluginObj.pan_to_selected = False
        pluginObj.plot(future2, plotObj)

        select_mode = p.select_mode.lower()
        manualSelect = False
        if select_mode != 'manual':
            if (not p.info.has_key('num_preferred') or 
                (p.info['num_preferred'] == 0) or 
                (len(p.starlist) == 0)):
                msg = msg_semiauto_failure
                self.fv.play_soundfile(snd_auto_failure, priority=19)
                manualSelect = True
            else:
                p.result = 'ok'
                p.selected = [ p.starlist[0] ]
                msg = msg_auto_success

        else:
            msg = msg_auto_manual
            manualSelect = True

        self.fv.ds.raise_tab('Dialogs')
        self.fv.ds.raise_tab('DSS: AgAutoSelect')
        pluginObj.set_message(msg)

        if not manualSelect:
            future2.resolve(0)
            return

        self.fv.update_pending(timeout=0.10)
        self.fv.play_soundfile(snd_auto_select_manual, priority=20)

    def _ag_auto_select_cont3(self, future2, future):
        self.logger.debug("continuation 3 resumed...")
        p = future.get_data()
        try:
            if p.result == 'ok':
                star = p.selected[0]
                star_ra_deg = star['ra_deg']
                star_dec_deg = star['dec_deg']
                star_mag = star['mag']
                star_name = star['name']

                # Return coords of picked star & magnitude
                p.star_ra = radec.raDegToString(star_ra_deg)
                p.star_dec = radec.decDegToString(star_dec_deg)
                p.star_mag = star_mag
                p.star_name = star_name

                sep_ra, sep_dec = wcs.get_RaDecOffsets(star_ra_deg, star_dec_deg,
                                                       p.probe_ra_deg, p.probe_dec_deg)
                p.ra_off = radec.offsetRaDegToString(sep_ra)
                p.dec_off = radec.decDegToString(sep_dec)

                # Calculate exposure region based on foci
                expregion = ag_config.calc_expregion(p.f_select)
                p.exp_x1 = expregion[0]
                p.exp_y1 = expregion[1]
                p.exp_x2 = expregion[2]
                p.exp_y2 = expregion[3]

                # Calculate exposure time based on magnitude
                p.exp_time = ag_config.calc_exposure(p.f_select, star_mag)

                agcodes = ag_config.calc_agcodes(p.f_select, p.ag_area)
                p.ag_x1 = agcodes[0]
                p.ag_y1 = agcodes[1]
                p.ag_x2 = agcodes[2]
                p.ag_y2 = agcodes[3]

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        # These won't pass back over remoteObjects
        p.setvals(starlist=None, selected=None, vignette_map=None,
                  image=None)
        
        self.logger.debug("ag_auto_select cb terminating: res=%s" % (str(p)))
        future.resolve(0)


    def sh_auto_select(self, tag, future,
                       motor=None, equinox=None, ra=None, dec=None,
                       select_mode=None, region=None):

        chname = 'DSS'
        if not self.fv.has_channel(chname):
            self.fv.add_channel(chname)
        chinfo = self.fv.get_channelInfo(chname)

        cat_fov = region          # is this in degrees?
        magic_constant = 2.5
        ra_deg = radec.funkyHMStoDeg(ra)
        dec_deg = radec.funkyDMStoDeg(dec)

        #f_select = 'P_OPT'
        #outer_fov = self.calculate_outer_fov(f_select)
        # NOTE: this seems to be closest to the fov shown in SOSS ShAutoSelect
        outer_fov = ag_config.calc_instrument_fov('SPCAM')

        dss_fov = region * magic_constant
        # Create blank image to load and calculate WCS for plotting
        px_scale = 0.00488281
        image = dp.create_blank_image(ra_deg, dec_deg, dss_fov,
                                      px_scale, 0.0,
                                      cdbase=[-1, 1],
                                      logger=self.logger)

        # Load image into DSS channel
        fitsname = 'SH_DUMMY'
        metadata = image.get_metadata()
        self.update_image(fitsname, chinfo.name, image, metadata)
        
        p = future.get_data()
        p.setvals(equinox=equinox, chname=chname, select_mode=select_mode,
                  image=image)
        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._sh_auto_select_cont2, future)
        
        # Open up the UI
        if not chinfo.opmon.is_active('AgAutoSelect'):
            chinfo.opmon.start_plugin_future('DSS', 'AgAutoSelect',
                                             future2, alreadyOpenOk=True)
        else:
            self.fv.ds.raise_tab('DSS')

        # calculate center pixel for ra/dec
        p.ctr_x, p.ctr_y = image.radectopix(ra_deg, dec_deg)

        # calculate radius of parameter degree radius fov
        radii = []
        r = 0.5
        while r <= region:
            cat_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, r)
            radii.append(cat_radius)
            r += 0.5
        p.cat_radii = radii

        # calculate radius of probe outer movable area fov
        p.outer_radius = image.calc_radius_xy(p.ctr_x, p.ctr_y, outer_fov)
        self.logger.info("Probe outer movable area radius is %d pixels." % (
            p.outer_radius))

        def query_catalogs(p):
            try:
                # Get preferred guide star catalog for SH
                catname = self.settings.get('SH_catalog', 'sh@subaru')
                starcat = self.catalogs.getCatalogServer(catname)

                radius = cat_fov * 60.0
                # Query catalog
                # query_result = starcat.search(
                #     ra_deg=ra_deg, dec_deg=dec_deg, equinox=2000.0,
                #     fov_deg=cat_fov, upper_mag=13.0)
                query_result = starcat.search(
                    ra=ra_deg, dec=dec_deg, equinox=2000.0,
                    r1=0.0, r2=radius, m1=0.0, m2=13.0)

                info, starlist = starcat.process_result(query_result)
                p.info = info
                self.logger.debug("info=%s" % (str(info)))
                p.starlist = starlist
                if len(starlist) > 0:
                    p.selected = [ starlist[0] ]
                else:
                    p.selected = []

            except Exception, e:
                errmsg = "Error querying star catalog: %s" % (str(e))
                self.logger.error(errmsg)
                self.fv.show_error(errmsg)
                p.setvals(info={}, starlist=[], selected=[], error=errmsg)
                raise VGWError(errmsg)

        # Offload this network task to a non-gui thread
        f_cat = Future.Future()
        f_cat.freeze(query_catalogs, p)

        chinfo = self.fv.get_channelInfo(chname)
        # Clear old data from canvas
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        pluginObj.reset()
        
        chinfo.fitsimage.onscreen_message("Querying catalog db...",
                                          delay=1.0)
        self.fv.showStatus("Querying catalog for objects ...")
        f_rest = Future.Future()
        f_rest.freeze(self.fv.gui_do, self._sh_auto_select_cont1, future2,
                      future)
        f_cat.add_callback('resolved', f_rest.thaw)
        self.fv.nongui_do_future(f_cat)

    def _sh_auto_select_cont1(self, future2, future):
        self.logger.debug("continuation 1 resumed...")
        p = future.get_data()

        if p.starlist is None:
            self.fv.show_error("No catalog data returned!")
            future.resolve(-1)
            return

        self.fv.showStatus("catalog returned %d results" % (
            len(p.starlist)))
        
        # Plot all the data
        plotObj = SHfov(self.logger, p.image, p)
        chinfo = self.fv.get_channelInfo(p.chname)
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        # TEMP: fix
        pluginObj.limit_stars_to_area = False
        pluginObj.pan_to_selected = False
        pluginObj.plot(future2, plotObj)

        select_mode = p.select_mode.lower()
        manualSelect = False
        if select_mode != 'manual':
            if not p.info.has_key('num_preferred') or \
               (p.info['num_preferred'] == 0):
                msg = msg_semiauto_failure
                manualSelect = True
                self.fv.play_soundfile(snd_auto_failure, priority=19)
            else:
                p.result = 'ok'
                p.selected = [ p.starlist[0] ]
                msg = msg_auto_success
        else:
            msg = msg_auto_manual
            manualSelect = True

        pluginObj.set_message(msg)

        if not manualSelect:
            future2.resolve(0)
            return

        self.fv.play_soundfile(snd_auto_select_manual, priority=20)

    def _sh_auto_select_cont2(self, future2, future):
        self.logger.debug("continuation 2 resumed...")
        p = future.get_data()
        try:
            if p.result == 'ok':
                star = p.selected[0]
                star_ra_deg = star['ra_deg']
                star_dec_deg = star['dec_deg']
                star_mag = star['mag']
                star_name = star['name']

                # Return coords of picked star & magnitude
                p.ra = radec.raDegToString(star_ra_deg)
                p.dec = radec.decDegToString(star_dec_deg)
                p.mag = star_mag
                p.name = star_name

                # Calculate exposure time based on magnitude
                p.exp_ag = ag_config.calc_exposure_sh('ag', star_mag)
                p.exp_sh = ag_config.calc_exposure_sh('sh', star_mag)

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        # These won't pass back over remoteObjects
        p.setvals(starlist=None, selected=None, vignette_map=None,
                  image=None)

        self.logger.debug("sh_auto_select cb terminating: res=%s" % (str(p)))
        future.resolve(0)


    def hsc_ag_auto_select(self, tag, future,
                           motor=None, equinox=None, ra=None, dec=None,
                           ag_pa=None, select_mode=None,
                           dss_mode=None, ag_area=None, instrument_name=None,
                           limitmag=None, goodmag=None, fov_pattern=None,
                           ut1_utc=None):

        chname = 'DSS'
        if not self.fv.has_channel(chname):
            self.fv.add_channel(chname)
        chinfo = self.fv.get_channelInfo(chname)

        ra_deg = radec.funkyHMStoDeg(ra)
        dec_deg = radec.funkyDMStoDeg(dec)

        # Calculate fov from instrument
        #inst_fov = 2.5
        inst_fov = ag_config.calc_instrument_fov(instrument_name)
        
        # For DSS, ra and dec are specified in traditional format
        ra_txt = radec.raDegToString(ra_deg, format='%02d:%02d:%06.3f')
        dec_txt = radec.decDegToString(dec_deg,
                                       format='%s%02d:%02d:%05.2f')

        if not dss_mode:
            dss_mode = 'off'
        else:
            dss_mode = dss_mode.lower()
            if dss_mode == 'on':
                dss_mode = 'dss1'
                
        # DSS FOV (deg)
        dss_fov_deg = inst_fov * 1.66666
        # Save some information for the continuation
        p = future.get_data()
        p.setvals(image=None, ra_deg=ra_deg, dec_deg=dec_deg,
                  inst_fov=inst_fov, ag_pa=ag_pa, select_mode=select_mode,
                  dss_mode=dss_mode, limitmag=limitmag, goodmag=goodmag,
                  ut1_utc=ut1_utc, chname=chname, equinox=equinox, ag_area=ag_area)

        future2 = Future.Future(data=p)
        future2.add_callback('resolved', self._hsc_ag_auto_select_cont3,
                             future)
        
        # Open up the UI
        if not chinfo.opmon.is_active('AgAutoSelect'):
            chinfo.opmon.start_plugin_future(chname, 'AgAutoSelect',
                                             future2, alreadyOpenOk=True)
        else:
            self.fv.ds.raise_tab('DSS')

        def get_dss_image(p):
            pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
            pluginObj = pluginInfo.obj

            # Assume square image?
            wd_deg = dss_fov_deg
            ht_deg = dss_fov_deg
        
            # width and height are specified in arcmin
            sgn, deg, mn, sec = radec.degToDms(wd_deg)
            wd = deg*60.0 + float(mn) + sec/60.0
            sgn, deg, mn, sec = radec.degToDms(ht_deg)
            ht = deg*60.0 + float(mn) + sec/60.0

            # Get the DSS image into the window
            try:
                if dss_mode != 'off':
                    params = dict(ra=ra_txt, dec=dec_txt, width=wd, height=ht)

                    # Get image server name
                    default_server = dss_mode + self.dss_server_sfx
                    server = self.settings.get('dss_server', default_server)

                    # Query the server and download file
                    fitspath = pluginObj.get_sky_image(server, params)
                    image = self.fv.load_file(fitspath, chname=chinfo.name)
                else:
                    # make blank image
                    #px_scale = 0.00488281
                    px_scale = 0.000280178318866
                    image = dp.create_blank_image(ra_deg, dec_deg,
                                                  dss_fov_deg,
                                                  px_scale, 0.0,
                                                  cdbase=[-1, 1],
                                                  logger=self.logger)
                    self.fv.gui_call(chinfo.fitsimage.set_image, image)

                p.image = image

            except Exception as e:
                errmsg = "Error querying dss server: %s" % (str(e))
                p.image = None
                self.logger.error(errmsg)
                self.fv.show_error(errmsg)

        # Offload this network task to a non-gui thread
        f_dss = Future.Future()
        f_dss.freeze(get_dss_image, p)

        # Clear old data from canvas
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        pluginObj.reset()
        
        chinfo.fitsimage.onscreen_message("Querying image db...",
                                          delay=1.0)
        self.fv.showStatus("Querying sky image database ...")
        f_rest = Future.Future()
        f_rest.freeze(self.fv.gui_do, self._hsc_ag_auto_select_cont1, future2,
                      future)
        f_dss.add_callback('resolved', f_rest.thaw)
        self.fv.nongui_do_future(f_dss)

    def _hsc_ag_auto_select_cont1(self, future2, future):
        self.logger.debug("continuation 1 resumed...")
        p = future.get_data()
        if p.image is None:
            # TODO: pop up an error message
            self.fv.show_error("No DSS image returned!")
            future.resolve(-1)
            return
        self.fv.showStatus("Got DSS image.")

        image = p.image
        # Now that we have an image, we can do some WCS calculations
        
        # Kawanomoto's helper object needs these times
        ut1 = p.ut1_utc
        sec = time.time()
        lst_deg = (astro_wcs.calcLST_sec(sec, ut1) / 3600.0) * 15.0
        mjd = astro_wcs.calcMJD(sec, ut1)

        # Get the ra/dec positions of the guiding CCD corners
        coords = self.hscguideinfo.get_ccdpos(p.ra_deg, p.dec_deg, p.ag_pa,
                                              lst_deg, mjd)

        if p.ag_area == 'DITHER':
            agarea_coords = self.hscguideinfo.get_dither_ccdpos(p.ra_deg, p.dec_deg, p.ag_pa, lst_deg, mjd)
        else:  # ag_area == SINGLE
            agarea_coords = self.hscguideinfo.get_vignette_ccdpos(p.ra_deg, p.dec_deg,
                                                            p.ag_pa, lst_deg, mjd)

        # TEMP: there is a bad channel in CCD 2_18
        # so remove it from consideration for now
        ## bad_ccd_idx = 2
        ## coords.pop(bad_ccd_idx)
        ## agarea_coords.pop(bad_ccd_idx)

        # For each CCD, get the coordinates of the corners, accounting
        # for distortion, so we can draw them on the star field
        agarea_polygons = []
        agarea_pixel_polygons = []
 
        for i in range(len(agarea_coords)):
            corner = agarea_coords[i]

            # Find the middle of this CCD on the sky
            ra0, dec0 = corner[0]
            x0, y0 = image.radectopix(ra0, dec0)
            ra1, dec1 = corner[1]
            x1, y1 = image.radectopix(ra1, dec1)
            ra2, dec2 = corner[2]
            x2, y2 = image.radectopix(ra2, dec2)
            ra3, dec3 = corner[3]
            x3, y3 = image.radectopix(ra3, dec3)

            points = [(ra0, dec0), (ra1, dec1), (ra2, dec2), (ra3, dec3)]
            agarea_polygons.append(points)

            points = [(x0, y0), (x1, y1), (x2, y2), (x3, y3)]
            agarea_pixel_polygons.append(points)

        queries = []
        polygons = []
        circles = []
        for i in range(len(coords)):
            corner = coords[i]

            # Find the middle of this CCD on the sky
            ra0, dec0 = corner[0]
            x0, y0 = image.radectopix(ra0, dec0)
            ra1, dec1 = corner[1]
            x1, y1 = image.radectopix(ra1, dec1)
            ra2, dec2 = corner[2]
            x2, y2 = image.radectopix(ra2, dec2)
            ra3, dec3 = corner[3]
            x3, y3 = image.radectopix(ra3, dec3)

            points = [(x0, y0), (x1, y1), (x2, y2), (x3, y3)]
            self.logger.debug("Points for ccd %d: %s" % (i, str(points)))
            polygons.append(points)

            x, y = (x0 + x2)/2, (y0 + y2)/2
            ra, dec = image.pixtoradec(x, y)
            self.logger.debug("Center for ccd %d: x,y=%f,%f ra,dec=%f,%f" % (
                i, x, y, ra, dec))

            # Find the radius for a circular search that will give us all
            # the guide stars in the rectangular CCD
            radius_pix = math.sqrt(math.fabs(x2 - x)**2 + math.fabs(y2 - y)**2)
            radius_deg = wcs.deltaStarsRaDecDeg(ra, dec, ra2, dec2)
            radius = radius_deg * 60.0
            circles.append((x, y, radius_pix))
            queries.append((ra, dec, radius))

        p.setvals(polygons=polygons, circles=circles,
                  starlist=None, agarea_polygons=agarea_polygons, 
                  agarea_pixel_polygons=agarea_pixel_polygons,
                  queries=queries)

        #agarea_polygons=agarea_polygonsw
        def query_catalogs(queries, p):
            # TODO: make these database searches concurrent
            all_stars = []
            try:
                # Get preferred guide star catalog for HSC
                catname = self.settings.get('HSC_catalog', 'hscag@subaru')
                starcat = self.catalogs.getCatalogServer(catname)

                # Query catalog
                for (ra, dec, radius) in queries:
                    self.logger.debug("Querying star catalog (ccd %d): ra=%f dec=%f r=%f" % (
                        i, ra, dec, radius_deg))
                    starlist, info = starcat.search(
                        ra=str(ra), dec=str(dec), r1=str(0.0), r2=str(radius),
                        catalog='', m2=str(p.limitmag), m1=str(p.goodmag))

                    all_stars.extend(starlist)
                    p.info = info
                    self.logger.debug("info=%s" % (str(info)))

            except Exception, e:
                errmsg = "Error querying star catalog: %s" % (str(e))
                p.setvals(result='error', errmsg=errmsg)
                self.logger.error(errmsg)
                self.fv.show_error(errmsg)
                p.setvals(info={}, starlist=[], selected=[], error=errmsg,
                          image=None, queries=None, circles=None, polygons=None,
                          agarea_pixel_polygons=None, agarea_polygons=None)
                raise VGWError(errmsg)

            p.starlist = all_stars
            if len(all_stars) > 0:
                p.selected = [ all_stars[0] ]
            else:
                p.selected = []

        # Offload this network task to a non-gui thread
        f_cat = Future.Future()
        f_cat.freeze(query_catalogs, queries, p)

        chinfo = self.fv.get_channelInfo(p.chname)
        chinfo.fitsimage.onscreen_message("Querying catalog db...",
                                          delay=1.0)
        self.fv.showStatus("Querying catalog for objects ...")
        f_rest = Future.Future()
        f_rest.freeze(self.fv.gui_do, self._hsc_ag_auto_select_cont2, future2,
                      future)
        f_cat.add_callback('resolved', f_rest.thaw)
        self.fv.nongui_do_future(f_cat)

    def _hsc_ag_auto_select_cont2(self, future2, future):
        self.logger.debug("continuation 2 resumed...")
        p = future.get_data()

        if p.starlist is None:
            self.fv.show_error("No catalog data returned!")

            future.resolve(-1)
            return

        self.fv.showStatus("catalog returned %d results" % (
            len(p.starlist)))
        
        # Plot all the data
        plotObj = HSCfov(self.logger, p.image, p)

        chinfo = self.fv.get_channelInfo(p.chname)
        pluginInfo = chinfo.opmon.getPluginInfo('AgAutoSelect')
        pluginObj = pluginInfo.obj
        # TEMP: fix
        pluginObj.limit_stars_to_area = True
        pluginObj.pan_to_selected = True
        try:
            pluginObj.plot(future2, plotObj)

        except Exception, e:
            errmsg = "Error filtering stars: %s" % (str(e))
            p.setvals(result='error', errmsg=errmsg)
            self.logger.error(errmsg)
            self.fv.show_error(errmsg)
            p.setvals(info={}, starlist=[], selected=[], error=errmsg,
                      image=None, queries=None, circles=None, polygons=None,
                      agarea_pixel_polygons=None, agarea_polygons=None)
            future.resolve(-1)
            return

        # scale image to 100% for easy viewing
        chinfo.fitsimage.scale_to(1.0, 1.0, no_reset=True)

        select_mode = p.select_mode.lower()
        manualSelect = False
        if select_mode != 'manual':
            if not p.info.has_key('num_preferred') or \
               (p.info['num_preferred'] == 0):
                msg = msg_semiauto_failure
                self.fv.play_soundfile(snd_auto_failure, priority=19)
                manualSelect = True
            else:
                p.result = 'ok'
                p.selected = [ p.starlist[0] ]
                msg = msg_auto_success

        else:
            msg = msg_auto_manual
            manualSelect = True

        pluginObj.set_message(msg)

        if not manualSelect:
            future2.resolve(0)
            return

        self.fv.update_pending(timeout=0.10)
        self.fv.play_soundfile(snd_auto_select_manual, priority=20)

    def _hsc_ag_auto_select_cont3(self, future2, future):
        self.logger.debug("continuation 3 resumed...")
        p = future.get_data()
        try:
            if p.result == 'ok':
                star = p.selected[0]
                self.logger.debug("selected star is %s" % str(star))
                star_ra_deg = star['ra_deg']
                star_dec_deg = star['dec_deg']
                star_mag = star['mag']
                star_name = star['name']

                # Return coords of picked star & magnitude
                p.star_ra = radec.raDegToString(star_ra_deg)
                p.star_dec = radec.decDegToString(star_dec_deg)
                p.star_mag = star_mag
                p.star_name = star_name

                # Calculate exposure time based on magnitude
                #p.exp_time = ag_config.calc_exposure('popt2', star_mag)
                p.exp_time = 0

        except Exception, e:
            p.setvals(result='error', errmsg=str(e))

        # These won't pass back over remoteObjects
        p.setvals(starlist=None, selected=None, queries=None,
                  image=None, circles=None, polygons=None,
                  agarea_pixel_polygons=None, agarea_polygons=None)
        
        self.logger.debug("hsc_ag_auto_select cb terminating: res=%s" % (str(p)))
        future.resolve(0)

    #############################################################
    #    Helper methods
    #############################################################

    def update_image(self, name, chname, image, metadata):
        """Called by the guider subsystem to place a guide image into the
        viewer.
        """
        self.fv.assert_gui_thread()

        if not self.fv.has_channel(chname):
            self.fv.add_channel(chname)

        chinfo = self.fv.get_channelInfo(chname)
        fitsimage = chinfo.fitsimage
        
        # Actually update the image
        self.fv.add_image(name, image, chname=chname)

        readout = chinfo.readout
        if readout:
            recv_time = time.time()
            frac_part = str(recv_time - int(recv_time)).split('.')[1][:3]
            time_str = "%s.%s" % (time.strftime("%H:%M:%S"), frac_part)
            readout.set_text(time_str)
            
        calctag = "%s-calc" % (chname)
        try:
            obj = fitsimage.getObjectByTag(calctag)

        except KeyError:
            obj = None

        # Check if we are guiding, and if so, show the calculation region
        # as an annotated rectangle
        if metadata.has_key('guiding'):
            if metadata['guiding']:
                x1, y1, x2, y2 = metadata['region']
                self.logger.debug("Guiding is on region is (%d,%d) (%d,%d)" % (
                    x1, y1, x2, y2))
                if obj is None:
                    fitsimage.add(CanvasTypes.CompoundObject(
                        CanvasTypes.Rectangle(x1, y1, x2, y2,
                                              color=self.colorcalc),
                        CanvasTypes.Text(x1, y2+4, "Calc Region",
                                         color=self.colorcalc),
                        ), tag=calctag, redraw=False)
                else:
                    bbox = obj.objects[0]
                    bbox.x1, bbox.y1, bbox.x2, bbox.y2 = x1, y1, x2, y2
            else:
                # Not guiding, clear the region, if any
                if obj is not None:
                    fitsimage.deleteObject(obj)

        
    def display_fitsbuf(self, fitsname, chname, data, width, height, na_type,
                        header, metadata):
        """Display a FITS image buffer.  Parameters:
        _fitsfile_: name of the file or title for the title bar
        _data_: ascii encoded numpy containing image data
        _width_, _height_: image dimensions in pixels
        _na_type_: numpy data type (currently ignored)
        _header_: fits file header as a dictionary
        _metadata_: metadata about image to attach to image
        """

        # Decode binary data
        data = ro.binary_decode(data)
        self.logger.debug("Received data: len=%d width=%d height=%d" % (
            len(data), width, height))
        self.logger.debug("metadata=%s header=%s" % (metadata, header))

        # Temporarily coerce numpy type  TODO: fix
        try:
            na_type = numpy.float32
            data = numpy.fromstring(data, dtype=na_type)
            data.byteswap(True)
            data = data.reshape((height, width))
            #print data

        except Exception, e:
            # Some kind of error decoding the value
            self.logger.error("Error creating image data for '%s': %s" % (
                fitsname, str(e)))
            return ro.ERROR

        # Create image container
        image = AstroImage.AstroImage(data, metadata=metadata,
                                      logger=self.logger)
        image.set(name=fitsname)
        image.update_keywords(header)
        
        # Enqueue image to display datasrc
        #self.fv.gui_do(self.fv.add_image, fitsname, image,
        #                    chname=chname)
        self.fv.gui_do(self.update_image, fitsname, chname, image,
                       metadata)

        return ro.OK

    def copy_data(self, chname1, chname2):
        try:
            chinfo1 = self.fv.get_channelInfo(chname1)
            chinfo2 = self.fv.get_channelInfo(chname2)

            image = chinfo1.fitsimage.get_image()
            chinfo2.fitsimage.set_image(image)
            return image

        except Exception, e:
            errmsg = "Failed to load %s image into %s: is there an image present?" % (
                chname1, chname2)
            self.logger.error(str(e))
            self.logger.error(errmsg)
            raise VGWError(errmsg)

    def map_back_to_ccd(self, p, image):
        self.logger.info("Region selection returned %s" % p)
        agh = Bunch.Bunch(image.get('agheader'))
        
        # Convert pixel coords on image back to CCD coords
        x1, y1 = self.pix2ccd(p.x1, p.y1, agh.binX,
                              agh.expRangeX, agh.expRangeY)
        x2, y2 = self.pix2ccd(p.x2, p.y2, agh.binX,
                              agh.expRangeX, agh.expRangeY)
        obj_x, obj_y = self.pix2ccd(p.obj_x, p.obj_y, agh.binX,
                              agh.expRangeX, agh.expRangeY)
        self.logger.info("CCD mapping is x1,y1=%d,%d x2,y2=%d,%d obj_x,obj_y=%f,%f" % (
            x1, y1, x2, y2, obj_x, obj_y))
        p.ccd_x1 = x1
        p.ccd_y1 = y1
        p.ccd_x2 = x2
        p.ccd_y2 = y2
        p.ccd_objx = obj_x
        p.ccd_objy = obj_y

    def pix2ccd(self, xi, yi, iBin, xoff, yoff):
        ## xo = xoff + ((xi-1) * iBin)
        ## yo = yoff + ((yi-1) * iBin)
        xo = xoff + (xi * iBin)
        yo = yoff + (yi * iBin)
        self.logger.info("xo=%f yo=%f xoff=%d yoff=%d iBin=%d xi=%f yi=%f" % (
            xo, yo, xoff, yoff, iBin, xi, yi))
        return (xo, yo)

    def ccd2pix(self, xi, yi, iBin, xoff, yoff):
        # xo = float(xi - xoff) / float(iBin)
        # yo = float(yi - yoff) / float(iBin)
        xo = (xi - xoff) // iBin
        yo = (yi - yoff) // iBin
        return (xo, yo)


    def _mark(self, chname, canvas, x, y, mode, mark, size, color):
        # mode is CLEAR | DRAW
        mode = mode.upper()
        if mode == 'CLEAR':
            objs = canvas.getObjectsByTagpfx("vgw_mark")
            canvas.deleteObjects(objs)
            
        elif mode == 'DRAW':
            color = color.lower()
            mark  = mark.lower()
            self.count += 1
            tag = ("vgw_mark%d" % self.count)

            # mark is POINT | CROSS | CIRCLE | SQUARE
            if mark == 'cross':
                tag = canvas.add(CanvasTypes.Point(x, y, size,
                                                    color=color),
                                 tag=tag)
            elif mark == 'point':
                tag = canvas.add(CanvasTypes.Circle(x, y, size,
                                                    color=color,
                                                    fill=True),
                                 tag=tag)
            elif mark == 'circle':
                tag = canvas.add(CanvasTypes.Circle(x, y, size,
                                                    color=color),
                                 tag=tag)
            elif mark == 'square':
                half = size #// 2
                x1, y1, x2, y2 = x-half, y-half, x+half, y+half
                tag = canvas.add(CanvasTypes.Rectangle(x1, y1, x2, y2,
                                                       color=color),
                                 tag=tag)

            # Only raise for a draw
            self.fv.ds.raise_tab(chname)

        
    def __str__(self):
        return 'vgw'
    

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
        
        return map(mm2pix, vlist)

    def draw(self, pluginObj):
        canvas = pluginObj.get_canvas()
        p = self.p
        thickness = self.ring_thickness
        color = self.colors

        # Draw AG probe outer movable area fov
        canvas.add(
            CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.outer_radius,
                               color=color.outer,
                               linestyle='dash', linewidth=thickness),
            redraw=False)

        ## # Draw AG probe inner movable area fov
        ## canvas.add(
        ##     CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.inner_radius,
        ##                        color=color.inner,
        ##                        linestyle='dash', linewidth=thickness),
        ##     redraw=False)

        # Draw AG probe position as a circle
        canvas.add(
            CanvasTypes.Circle(p.probe_x, p.probe_y, self.probe_radius,
                               color=color.probe,
                               linestyle='dash', linewidth=thickness),
            redraw=False)

        # Draw vignette map
        self.vig_obj = CanvasTypes.Polygon(self.vignette_map,
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
            CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.inst_radius,
                               color=color.inst,
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
        
        obj = CanvasTypes.CompoundObject()
        canvas.add(obj, redraw=False)

        obj.addObject(CanvasTypes.Polygon(self.fov_pts,
                                          color=color.inst,
                                          linestyle='dash',
                                          linewidth=thickness))
        obj.addObject(CanvasTypes.Polygon(self.vig_pts,
                                            color='blue',
                                            linestyle='solid',
                                            linewidth=thickness))
        obj.addObject(CanvasTypes.Text(self.c1x, self.c1y, "Chip1",
                                       color='white'))
        obj.addObject(CanvasTypes.Text(self.c2x, self.c2y, "Chip2",
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
        
        obj = CanvasTypes.CompoundObject()
        canvas.add(obj, redraw=False)

        # Draw instrument fov
        obj.addObject(
            CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.inst_radius,
                               color=color.inst,
                               linestyle='dash', linewidth=thickness))

        # Draw outer telescope foci fov
        obj.addObject(
            CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.outer_radius,
                               color=color.outer,
                               linestyle='dash', linewidth=thickness))

        # Draw probe position
        obj.addObject(
            CanvasTypes.Circle(p.probe_x, p.probe_y, self.probe_radius,
                               color=color.probe,
                               linestyle='dash', linewidth=thickness))
        ## # Draw AG probe inner movable area fov
        ## obj.addObject(
        ##     CanvasTypes.Circle(p.ctr_x, p.ctr_y, self.inner_radius,
        ##                        color=color.inner,
        ##                        linestyle='dash', linewidth=thickness))

        # Draw probe movable area
        self.prb_area = CanvasTypes.Polygon(self.fov_pts,
                                            color=color.probe,
                                            linestyle='dash',
                                            linewidth=thickness)
        obj.addObject(self.prb_area)

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
            gons.append(CanvasTypes.Polygon(points,
                                            color=color.inst,
                                            linestyle='dash',
                                            linewidth=thickness))
        self.ccds_obj = CanvasTypes.CompoundObject(*gons)
        canvas.add(self.ccds_obj, redraw=False)

        pixgons = []
        for points in p.agarea_pixel_polygons:
            # Add internal dither image polygon
            self.logger.info("Plotting dithering polygon %s" % str(points))
            pixgons.append(CanvasTypes.Polygon(points,
                                            color=color.vignette,
                                            linestyle='dash',
                                            linewidth=1))
        self.dith_obj = CanvasTypes.CompoundObject(*pixgons)
        canvas.add(self.dith_obj, redraw=False)

        # for (x, y, r) in p.circles:
        #     canvas.add(
        #         CanvasTypes.Circle(x, y, r,
        #                            color=color.probe,
        #                            linestyle='solid', linewidth=1),
        #         redraw=False)

    def filter_results(self, starlist):
        # return self.pluginObj.filter_results(starlist, self.dith_obj)
        #all_stars = map(lambda stars: self.pluginObj.filter_results(stars,
        #                   self.dith_obj), starlist)
        starlist = self.pluginObj.filter_results(starlist, self.dith_obj)

        self.logger.info('all STARLIST=%s' %str(starlist))

        p = self.p
        #return self.hsc_filter_candidates(p, p.queries, all_stars)
        return self.hsc_filter_candidates2(p, p.queries, starlist)

    def _set_priority(self, updated_stars):

        for num, star in enumerate(updated_stars):
            num += 1
            star['priority'] = num

    def hsc_filter_candidates2(self, p, queries, all_stars):
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

        all_stars = sorted(all_stars, key=itemgetter('preference')) 
        self. _set_priority(all_stars)

        # for n, s in zip(xrange(20), all_stars):
        #     self.logger.debug('name=%s, mag=%s flag=%s pref=%s' %(s['name'], str(s['mag']), str(s['flag']), str(s['preference'])))
        return all_stars

    def hsc_filter_candidates(self, p, queries, all_stars):
        # interesting items:
        # p.ag_pa, p.limitmag, p.goodmag
        # p.polygons, p.circles
        # queries: ((ra, dec, radius), ... ) 1 for each ccd
        # all_stars: query result for each circle

        fabs = math.fabs

        goodmag = p['goodmag']
        bad_star = 10
        mag_weight = 0.6 # mag's weight. 60% of flag weight 

        for star in all_stars:
            flag = star['flag']
            if flag < 0: # negative flag
                flag = fabs(flag-bad_star)
            calc_mag = (fabs(star['mag']-goodmag)) * mag_weight
            star['preference'] = flag + calc_mag     
            #print '%s flag=%s mag=%s calcmag=%s pref=%s' %(star['name'], flag, star['mag'], calc_mag, pref)   
            #star['preference'] = pref

        all_stars = sorted(all_stars, key=itemgetter('preference')) 
        self. _set_priority(all_stars)

        return all_stars

        # updated_stars = []

        # self.logger.info('Ps=%s' %str(p))

        # for query, stars in zip(queries, all_stars):
        #     self.logger.info('query=%s' %str(query))
        #     center_rarad = math.radians(query[0])
        #     center_decrad = math.radians(query[1])

        #     for star in stars:
        #         star_rarad = math.radians(star['ra_deg'])
        #         star_decrad = math.radians(star['dec_deg']) 
        #         delta = radec.delta_stars(center_rarad, center_decrad,
        #                                     star_rarad, star_decrad) 
        #         star['preference'] = delta
        #     updated_stars.extend(stars)
        # updated_stars = sorted(updated_stars, key=itemgetter('preference')) 
        # self._set_priority(updated_stars)

        # p.num_preferred = min(15, len(updated_stars))
        # return updated_stars
      
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
            CanvasTypes.Circle(p.ctr_x, p.ctr_y, p.outer_radius,
                               color=color.outer,
                               linestyle='dash', linewidth=thickness),
            redraw=False)

        # Draw coincentric catalog radii every 0.5 deg
        for cat_radius in p.cat_radii:
            canvas.add(
                CanvasTypes.Circle(p.ctr_x, p.ctr_y, cat_radius,
                                   color='white',
                                   linestyle='dash', linewidth=thickness),
                redraw=False)


#END
