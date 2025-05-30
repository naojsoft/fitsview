#
# PFS_AG.py -- PFS AG routines
#
# E. Jeschke
#
# Q for Hiro:
#  ) Does the "guide objects" table contain all GAIA stars found within the
#    boundary region of all CCDs?  Or just within a specific mag range?
#    Will there be only objects in this table that are detected?
#    > The latter. That table is supposed to contain all the GAIA stars
#    > specified by a particular "PFS design" in the operational database.
#    > They were filtered by their magnitude and others.
#  ) What does "spot_id" refer to within the detected objects table?
#    > The identified object table ties the "detected" objects to the "guide"
#    > catalog objects. The detected objects table contains "spot id"s,
#    > which are sequential numbers assigned by the AG camera controller.
#  1) what does "guide object id" column reference in "identified objects"
#       table
#     (it doesn't seem to be the same as "source id" in the "guide objects"
#       table--is it the row number?)
#     > yes
#  ) centroid_x, centroid_y in the detected_objects table seems to match the
#     centers of the objects (which makes sense). However, "guide_object_xdet"
#     and "guide_object_ydet" seem to be offset from the centers--what do these
#     values indicate?
#     > They are expected positions of the matched guide objects - if the
#     > pointing and detected-guide object mapping are perfect, they should
#     > match with centroid xy's.
#  2) what does the "flags" column in "detected objects" table mean?
#     > This may change, but currently, the flag indicates whether the
#     > detected spot is on the (nominally) outside of the focus side (0)
#     > or the inside of the focus side (1). The half of the AG camera
#     > detectors is covered by glass and the cameras are positioned so that
#     > each halves are (nominally) +/- 150 um off-focus. This is to measure
#     > the focus of PFI (spot sizes will be equal on both sides when the
#     > instrument is in focus). We are thinking to add flag bits to indicate
#     > whether a spot is on the outside-/inside-focus boundary, too close to
#     > the edge of the detector, etc.
#  3) will it be the case that sometimes we do not have a guide image for
#      some camera (under normal working conditions)?
#     > Yes. We expect those HDUs will still contain real, close to raw,
#     > images out of the AG cameras. There may be some expected guide objects
#     > in the guide objects table, but nothing in the detected objects table,
#     > obviously.
#
# > Please note that the focal plane coordinates (detected_object_x,
# > detected_object_y, guide_object_x, guide_object_y in mm) currently
# > used are HSC's, not PFS'. The sign of Y values is flipped (the X-axis
# > in both conventions points toward the opt side; the Z-axis in PFS
# > convention points toward the sky, while the Z-axis in HSC convention
# > points toward the primary mirror; so the Y-axis of PFS points "up" or
# > toward "rear", while the Y-axis of HSC points "down" or toward "front"
# > of the telescope. I should flip the sign so that what observers will
# > see are PFS focal plane coordinates.

# stdlib
import sys
import os
import time
import threading
import tempfile

# 3rd party
import numpy as np
from fitsio import FITS
# pip install inotify
import inotify.adapters
import yaml

# ginga
from ginga import trcalc, cmap, imap
from ginga.gw import ColorBar, Widgets, Viewers
from ginga.AstroImage import AstroImage
from ginga.util.io.io_fits import FitsioFileHandler
from ginga.util.wcsmod.wcs_astropy import AstropyWCS
from ginga.misc import Bunch
from ginga import GingaPlugin

# g2cam
from g2cam.status.client import StatusClient

# local
from fitsview.util import pfswcs


class PFS_AG(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super().__init__(fv)

        # get PFSAG preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_PFS_AG')
        self.settings.add_defaults(plot_fov=True,
                                   plot_identified_stars=True,
                                   plot_detected_not_identified=True,
                                   plot_guide_stars=True,
                                   detected_color='yellow',
                                   guide_color='cyan',
                                   identified_color='orangered',
                                   plot_offsets=False,
                                   subtract_bias=False,
                                   subtract_background=False,
                                   #color_map='Greens',
                                   #intensity_map='ramp',
                                   channel_name='PFS_1K',
                                   rate_limit=5.0,
                                   data_directory='.',
                                   save_directory=tempfile.gettempdir(),
                                   auto_orient=False,
                                   in_gen2=True)
        self.settings.load(onError='silent')

        self.chname = self.settings.get('channel_name')
        self.wsname = 'PFS_AG_CAMS'
        self.wstype = 'grid'
        self._in_gen2 = self.settings.get('in_gen2', True)
        self.inspace = 'sub1'
        self.fov_chname = 'PFS_FOV'
        self.fov_inspace = 'sub2'

        self._wd = 300
        self._ht = 300
        # self.mag_max = 20.0
        # self.mag_min = 12.0
        self.current_file = None
        self.current_time = ''
        self.ev_quit = threading.Event()
        # self.dark = dict()
        # self.flat = dict()
        # self.cmap_names = list(cmap.get_names())
        # self.cmap = cmap.get_cmap(self.settings.get('color_map'))
        # self.imap_names = list(imap.get_names())
        # self.imap = imap.get_imap(self.settings.get('intensity_map'))
        # self.field_names = ['mag']
        # self.field = 'mag'
        self.pause_flag = False
        self.rate_limit = self.settings.get('rate_limit', 5.0)
        self.error_scale = 1.0
        self.save_dir = self.settings.get('save_directory', '/tmp')
        self.last = Bunch.Bunch(image_time=time.time(),
                                raw=None, time_raw=None,
                                processed=None, time_processed=None)
        self.mode = 'processed'
        self.guide_count = 0

        # hold tables of detected objs, guide objs and identified objs
        self.tbl_do = None
        self.tbl_go = None
        self.tbl_io = None
        self.img_dct = {}
        self._fov_coords = dict()

        self.viewer = dict()
        self.dc = fv.get_draw_classes()

        # controls how far from the center we plot the guide cam images
        # in the FOV viewer
        self.fov_radius = 1500
        # angles measured from X axis for cams 1-6
        self.fov_angles = [90, 150, 210, 270, 330, 30]
        # angles measured from Y axis for rotating images 1-6
        self.rot_angles = [0, 60, 120, 180, 240, 300]

        if not self._in_gen2:
            self.fv.add_callback('add-image', self.incoming_data_cb)

        self.sc = None
        configfile = os.path.join(os.environ['CONFHOME'], 'status', 'status.yml')
        with open(configfile, 'r') as in_f:
            buf = in_f.read()
        self.config_d = yaml.safe_load(buf)

        self.gui_up = False

    def build_gui(self, container):

        top = Widgets.VBox()
        top.set_border_width(2)
        top.set_spacing(2)

        ws = self.fv.error_wrap(self.fv.add_workspace, self.wsname,
                                self.wstype, inSpace=self.inspace)
        ws.extdata.chpfx = 'CAM'
        # add grid of PFS AG viewers

        for i in range(0, 2):
            for j in range(0, 3):
                idx = i * 3 + j

                chname = 'CAM{}'.format(idx + 1)
                prefs = self.fv.get_preferences()
                settings = prefs.create_category(f'channel_{chname}')
                # don't generate thumbs for these channels
                settings.set(numImages=1, genthumb=False,
                             raisenew=False, focus_indicator=True)
                settings.load(onError='silent')
                channel = self.fv.add_channel(chname, settings=settings,
                                              workspace=self.wsname)
                viewer = channel.viewer

                # add a label, in case we are in grid view
                lbl = self.dc.Text(20, 30, chname, color='blue2',
                                   fontsize=12, font='sans',
                                   bgcolor='floralwhite', bgalpha=0.8,
                                   coord='window')
                lbl.crdmap = viewer.get_coordmap('window')
                canvas = viewer.get_canvas()
                canvas.add(lbl, redraw=False)

        # create "big picture" FOV channel
        chname = self.fov_chname
        settings = prefs.create_category(f'channel_{chname}')
        settings.set(numImages=1, genthumb=False, raisenew=False)
        settings.load(onError='silent')
        channel = self.fv.add_channel(chname, settings=settings,
                                      workspace=self.fov_inspace)
        viewer = channel.viewer
        canvas = self.dc.DrawingCanvas()
        viewer.get_canvas().add(canvas, redraw=False)
        canvas.set_surface(viewer)
        canvas.ui_set_active(True)
        self.fov_canvas = canvas

        captions = [('Plot FOV', 'checkbutton'),
                     # "Method:", 'label', 'method', 'combobox'),
                    ('Plot identified stars', 'checkbutton',
                     'Plot offsets', 'checkbutton',
                     'Auto orient', 'checkbutton'),
                    ('Plot guide stars', 'checkbutton',
                     'Plot NI detected stars', 'checkbutton'),
                    ('Subtract bias', 'checkbutton'),
                    # ('Subtract dark', 'checkbutton', 'Darks:', 'llabel',
                    #  'dark_frames', 'entryset'),
                    ('Subtract background', 'checkbutton'),
                    # ('Divide flat', 'checkbutton', 'Flats:', 'llabel',
                    #  'flat_frames', 'entryset'),
                    ('Pause', 'togglebutton', 'Rate Limit', 'spinfloat',
                     'Save', 'button'),
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        tf = self.settings.get('plot_fov', False)
        b.plot_fov.set_state(tf)
        b.plot_fov.set_tooltip("Plot images/stars in PFS FOV")
        b.plot_fov.add_callback('activated', self.toggle_plot_fov_cb)

        tf = self.settings.get('plot_identified_stars', False)
        b.plot_identified_stars.set_state(tf)
        b.plot_identified_stars.set_tooltip("Plot the identified guide stars")
        b.plot_identified_stars.add_callback('activated', self.toggle_plot_identified_stars_cb)

        tf = self.settings.get('plot_offsets', False)
        b.plot_offsets.set_state(tf)
        b.plot_offsets.set_tooltip("Plot the offsets from the identified guide stars")
        b.plot_offsets.add_callback('activated', self.toggle_plot_offsets_cb)

        tf = self.settings.get('auto_orient', False)
        b.auto_orient.set_state(tf)
        b.auto_orient.set_tooltip("Rotate images to match FOV plot")
        b.auto_orient.add_callback('activated', self.auto_orient_cb)

        tf = self.settings.get('plot_guide_stars', False)
        b.plot_guide_stars.set_state(tf)
        b.plot_guide_stars.set_tooltip("Plot all guide stars")
        b.plot_guide_stars.add_callback('activated',
                                           self.toggle_plot_guide_stars_cb)

        tf = self.settings.get('plot_detected_not_identified', False)
        b.plot_ni_detected_stars.set_state(tf)
        b.plot_ni_detected_stars.set_tooltip("Plot the detected stars that were not identified")
        b.plot_ni_detected_stars.add_callback('activated',
                                              self.toggle_plot_ni_detected_stars_cb)

        tf = self.settings.get('subtract_bias', False)
        b.subtract_bias.set_state(tf)
        b.subtract_bias.add_callback('activated', self.subtract_bias_cb)
        b.subtract_bias.set_tooltip("Subtract bias calculated from overscan region")
        top.add_widget(w, stretch=0)

        tf = self.settings.get('subtract_background', False)
        b.subtract_background.set_state(tf)
        b.subtract_background.add_callback('activated', self.subtract_bg_cb)
        b.subtract_background.set_tooltip("Subtract background calculated from image")
        top.add_widget(w, stretch=0)

        # tf = self.settings.get('subtract_dark', False)
        # b.subtract_dark.set_state(tf)
        # b.subtract_dark.add_callback('activated', self.subtract_dark_cb)
        # b.subtract_dark.set_tooltip("Subtract dark from frame")
        # b.dark_frames.add_callback('activated', self.set_darks_cb)
        # b.dark_frames.set_tooltip("Enter a file that contains dark frames")
        # #b.subtract_dark.set_enabled(False)
        # top.add_widget(w, stretch=0)

        # tf = self.settings.get('divide_flat', False)
        # b.divide_flat.set_state(tf)
        # b.divide_flat.add_callback('activated', self.divide_flat_cb)
        # b.divide_flat.set_tooltip("Divide frame by flat")
        # #b.divide_flat.set_enabled(False)
        # b.flat_frames.add_callback('activated', self.set_flats_cb)
        # b.flat_frames.set_tooltip("Enter a file that contains flat frames")

        b.pause.set_state(self.pause_flag)
        b.pause.add_callback('activated', self.pause_cb)
        b.pause.set_tooltip("Pause updates")
        b.pause.set_state(self.pause_flag)
        b.pause.add_callback('activated', self.pause_cb)
        b.pause.set_tooltip("Pause updates")
        b.rate_limit.set_limits(0.0, 60.0, incr_value=1.0)
        b.rate_limit.set_value(self.rate_limit)
        b.rate_limit.add_callback('value-changed', self.set_rate_limit_cb)
        b.rate_limit.set_tooltip("Set the rate limit for PFS AG files")

        b.save.add_callback('activated', self.save_file_cb)
        b.save.set_tooltip("Save current FITS object")

        top.add_widget(w, stretch=0)

        hbox = Widgets.HBox()
        btn1 = Widgets.RadioButton("Raw")
        btn1.set_state(self.mode == 'raw')
        btn1.add_callback('activated', self.set_mode_cb, 'raw')
        btn1.set_tooltip("Show raw frames")
        hbox.add_widget(btn1)

        btn2 = Widgets.RadioButton("Processed", group=btn1)
        btn2.set_state(self.mode == 'processed')
        btn2.add_callback('activated', self.set_mode_cb, 'processed')
        btn2.set_tooltip("Show processed frames")
        hbox.add_widget(btn2)

        fr = Widgets.Frame("AG Files")
        fr.set_widget(hbox)
        top.add_widget(fr, stretch=0)

        hbox = Widgets.HBox()
        hbox.add_widget(Widgets.Label("Received at:"), stretch=0)
        self.w.last_time = Widgets.Label('')
        hbox.add_widget(self.w.last_time, stretch=1)
        top.add_widget(hbox, stretch=0)

        # add CAM pan buttons
        hbox = Widgets.HBox()
        for i in range(1, 7):
            w = Widgets.Button(f"CAM{i}")
            hbox.add_widget(w, stretch=1)
            w.add_callback('activated', self.pan_cam_cb, i-1)
            w.set_tooltip(f"Pan to CAM{i} in the PFS_FOV viewer")
        w = Widgets.Button("ALL")
        hbox.add_widget(w, stretch=1)
        w.add_callback('activated', self.pan_zoom_fit_cb)
        w.set_tooltip(f"See all images the PFS_FOV viewer")
        top.add_widget(hbox, stretch=0)

        # stretch spacer
        top.add_widget(Widgets.Label(''), stretch=1)

        captions = [("Scaling:", 'label', 'error_scaling', 'hslider',
                     'error_scale', 'entryset'),
                    # ('Set minmax', 'button', 'minf', 'entry', 'maxf', 'entry'),
                    # ('color_map', 'combobox', 'intensity_map', 'combobox',
                    #  'field', 'combobox')
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.error_scaling.set_limits(0, 50, incr_value=1)
        i = int(max(0.0, min(50.0, (self.error_scale - 1) / 0.25)))
        b.error_scaling.set_value(i)
        b.error_scaling.add_callback('value-changed', self.set_error_scale_cb)
        b.error_scaling.set_tooltip("Change the error scaling")
        b.error_scale.set_tooltip("Set the error scaling precisely")
        b.error_scale.add_callback('activated', self.set_error_scale2_cb)
        b.error_scale.set_text("{:.2f}".format(self.error_scale))

        # b.minf.set_text(str(self.mag_min))
        # b.minf.add_callback('activated', self.set_minmax_cb)
        # b.maxf.set_text(str(self.mag_max))
        # b.maxf.add_callback('activated', self.set_minmax_cb)
        # b.set_minmax.add_callback('activated', self.set_minmax_cb)
        # b.set_minmax.set_tooltip("Set min/max values for color bar")

        # combobox = b.color_map
        # for name in self.cmap_names:
        #     combobox.append_text(name)
        # cmap_name = 'rainbow3'
        # try:
        #     index = self.cmap_names.index(cmap_name)
        # except Exception:
        #     index = self.cmap_names.index('gray')
        # combobox.set_index(index)
        # combobox.add_callback('activated', self.set_cmap_cb)

        # combobox = b.intensity_map
        # for name in self.imap_names:
        #     combobox.append_text(name)
        # imap_name = 'ramp'
        # try:
        #     index = self.imap_names.index(imap_name)
        # except Exception:
        #     index = self.imap_names.index('ramp')
        # combobox.set_index(index)
        # combobox.add_callback('activated', self.set_imap_cb)

        # combobox = b.field
        # for name in self.field_names:
        #     combobox.append_text(name)
        # field_name = 'mag'
        # try:
        #     index = self.field_names.index(field_name)
        # except Exception:
        #     index = self.field_names.index('mag')
        # combobox.set_index(index)
        # combobox.add_callback('activated', self.set_field_cb)

        top.add_widget(w, stretch=0)

        # self.cbar = ColorBar.ColorBar(self.logger, link=True)
        # self.cbar.set_cmap(self.cmap)
        # self.cbar.set_imap(self.imap)
        # rgbmap = self.cbar.get_rgbmap()
        # rgbmap.add_callback('changed', self.replot_stars)
        # # hack to set font size of this color bar
        # self.cbar.cbar.fontsize = 8

        # cbar_w = self.cbar.get_widget()
        # cbar_w.resize(-1, 32)
        # top.add_widget(cbar_w, stretch=0)

        #top.add_widget(Widgets.Label(''), stretch=1)

        btns = Widgets.HBox()
        btns.set_border_width(4)
        btns.set_spacing(3)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        #btn = Widgets.Button("Help")
        #btn.add_callback('activated', lambda w: self.help())
        #btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)

        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)
        self.gui_up = True

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def start(self):
        self.ev_quit.clear()
        if self._in_gen2:
            self.fv.nongui_do(self.watch_loop, self.ev_quit)

        self.connect_status()

    def stop(self):
        self.ev_quit.set()

        # delete channels created by this plugin
        for idx in range(0, 6):
            chname = 'CAM{}'.format(idx + 1)
            self.fv.delete_channel(chname)
        self.fv.delete_channel(self.fov_chname)

        # delete the workspace
        try:
            ws = self.fv.ds.get_ws(self.wsname)
            self.fv.delete_workspace(ws)
        except KeyError:
            pass

        self.gui_up = False

    def focus_cb(self, viewer, tf):
        if not tf:
            return
        chname = self.fv.get_channel_name(viewer)
        self.fv.change_channel(chname, raisew=True)

    def incoming_data_cb(self, fv, chname, image, info):
        if chname != 'PFS':
            return
        # TEMP: until we remove PFS_AG_NEW
        if not self.gui_up:
            return

        #self.fv.nongui_do(self.process_image, image)

        path = image.get('path', None)
        if path is not None:
            self.incoming_file(path)

    def incoming_file(self, filepath):
        filedir, filename = os.path.split(filepath)
        # save filepath for last received one of each type
        if filename.startswith('_'):
            self.last.setvals(raw=filepath, time_raw=time.time())
            if self.mode == 'processed':
                return
        else:
            self.last.setvals(processed=filepath, time_processed=time.time())
            if self.mode == 'raw':
                return

        self.fv.nongui_do(self.process_file, filepath, set_1k=True)

    def quick_data_reduce(self, image, name):
        data = image.get_data()
        wd, ht = image.get_size()[:2]
        # overscans - 24 pixels at the beginning and end of each row (X)
        # and first 9 rows (Y).
        ovsc_x, ovsc_y = 24, 9

        if self.settings.get('subtract_bias', False):
            rows, columns = data.shape[:2]
            overscan = data[0:ovsc_y + 1, :]
            medians = np.nanmedian(overscan, axis=0).reshape((1, columns))
            subtrahend = np.repeat(medians, rows, axis=0)
            data = data - subtrahend

        if self.settings.get('subtract_background', False):
            med = np.nanmedian(data[ovsc_y:ht, ovsc_x:wd - ovsc_x])
            data = data - med

        # if self.settings.get('subtract_dark', False):
        #     if len(self.dark) == 0:
        #         self.logger.error("No darks are loaded")
        #     else:
        #         data = data - self.dark[name]

        # if self.settings.get('divide_flat', False):
        #     if len(self.dark) == 0:
        #         self.logger.error("No flats are loaded")
        #     else:
        #         med = np.median(self.flat[name])
        #         data = (data / self.flat[name]) * med

        new_img = AstroImage(data_np=data)
        new_img.update_keywords(image.get_header())
        new_img.set(name=image.get('name'), nothumb=True)
        #new_img.set(name='PFS' + str(time.time()), nothumb=True)
        return new_img

    # def process_image(self, image):
    #     path = image.get('path', None)
    #     if path is None:
    #         return

    #     self.process_file(path, set_1k=True)

    def update_grid(self, img_dct):
        self.fv.assert_gui_thread()

        self.img_dct = img_dct
        self.fv.assert_gui_thread()
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                if cam_id in img_dct:
                    image = img_dct[cam_id]
                    image.set(tag=cam_id)
                    canvas = viewer.get_canvas()
                    channel.add_image(image)
                else:
                    viewer.clear()

        self.auto_orient()

        # show when this file was received
        self.w.last_time.set_text(self.current_time)

        self.fv.update_pending()

        if self._in_gen2:
            # increment the guide count
            self.guide_count += 1
            stat_d = {'VGW.PFS.AG.COUNT': self.guide_count}
            self.fv.nongui_do(self.fv.call_global_plugin_method,
                              'Gen2Int', 'store', [stat_d], {})

        if self.settings.get('plot_fov', False):
            self.plot_fov()

        if self.tbl_go is not None:
            self.plot_stars()
        else:
            self.clear_stars(redraw=True)

    def set_1k(self, image):
        self.fv.assert_gui_thread()

        try:
            channel = self.fv.get_channel_on_demand(self.chname)
            channel.add_image(image)
        except Exception as e:
            self.logger.error(f"error setting 1K image: {e}", exc_info=True)

    def get_fov_xy(self, cam_num, x, y):
        base_x, base_y, sub_x, sub_y = self._fov_coords[cam_num]
        # NOTE: invert Y
        x, y = x - sub_x, ((sub_y * 2) - y) - sub_y
        rot_ang_deg = self.rot_angles[cam_num]
        x, y = trcalc.rotate_pt(x, y, rot_ang_deg)
        #pos_x, pos_y = np.asarray(base_x + x).astype(int), np.asarray(base_y + y).astype(int)
        pos_x, pos_y = base_x + x, base_y + y
        return (pos_x, pos_y)

    def plot_fov(self):
        channel = self.fv.get_channel_on_demand(self.fov_chname)
        viewer = channel.viewer
        canvas = self.fov_canvas
        canvas.delete_all_objects()

        for image in self.img_dct.values():
            if image is not None:
                data_np = image.get_data()
                dtype = data_np.dtype
                break

        self._fov_coords = dict()

        with viewer.suppress_redraw:
            dim = 4600
            ctr_x, ctr_y = dim * 0.5, dim * 0.5
            dst_arr = np.zeros((dim, dim, 2), dtype=dtype)
            for i, theta in enumerate(self.fov_angles):
                rot_ang_deg = self.rot_angles[i]
                theta = 90.0 + rot_ang_deg
                cam_num = i
                cam_id = "CAM{}".format(i + 1)
                if cam_id not in self.img_dct:
                    # guide camera MIA
                    continue
                image = self.img_dct[cam_id]
                data_np = image.get_data()
                # NOTE: invert Y side of image for correct orientation
                data_np = np.flipud(data_np)
                _wd, _ht = image.get_size()
                # add alpha layer
                mn, mx = trcalc.get_minmax_dtype(dtype)
                data_np = trcalc.add_alpha(data_np, alpha=mx)
                # rotate into correct orientation
                rot_data_np = trcalc.rotate(data_np, -rot_ang_deg)
                ht, wd = rot_data_np.shape[:2]
                off_x = self.fov_radius * np.cos(np.radians(theta))
                off_y = self.fov_radius * np.sin(np.radians(theta))
                x, y = off_x - wd * 0.5, off_y - ht * 0.5

                # drop image into position in final image
                pos_x, pos_y = (int(ctr_x + x), int(ctr_y + y))
                self.logger.debug(f"placing {cam_id} image at {pos_x},{pos_y}")
                dst_arr[pos_y:pos_y+ht, pos_x:pos_x+wd, :] += rot_data_np[:, :, :]
                # mark center of each camera
                base_x, base_y = ctr_x + off_x, ctr_y + off_y
                self._fov_coords[cam_num] = (base_x, base_y, _wd * 0.5, _ht * 0.5)
                canvas.add(self.dc.Point(base_x, base_y, radius=20,
                                         color='green',
                                         style='plus', linewidth=1),
                           redraw=False)

                # label each cam
                radius_txt = _ht * 0.5 + 10
                txt_x = radius_txt * np.cos(np.radians(theta)) + base_x
                txt_y = radius_txt * np.sin(np.radians(theta)) + base_y
                canvas.add(self.dc.Text(txt_x, txt_y, text=cam_id,
                                        rot_deg=rot_ang_deg,
                                        fillcolor='green'),
                           redraw=False)

            # mark center position
            canvas.add(self.dc.Point(ctr_x, ctr_y, radius=20, color='pink',
                                     style='plus', linewidth=1),
                       tag='center', redraw=False)
            rot_image = AstroImage(data_np=dst_arr, logger=self.logger)
            rot_image.set(name="PFS_guider_all", nothumb=True, path=None)
            viewer.set_image(rot_image)
            viewer.redraw(whence=0)

    def process_file(self, path, set_1k=False):
        self.fv.assert_nongui_thread()

        start_time = time.time()

        if self.current_file is not None:
            self.remove_file(self.current_file)
        self.current_file = path

        # filename beginning with an underscore is supposedly a raw file
        # with no WCS and no summary detection tables
        _dir, fname = os.path.split(path)
        fname, ext = os.path.splitext(fname)
        is_raw = fname.startswith('_')
        self.logger.info(f"processing '{fname}', raw={is_raw} ...")
        _t = self.last.time_raw if is_raw else self.last.time_processed
        if _t is not None:
            self.current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(_t))

        self.tbl_go = None
        self.tbl_do = None
        self.tbl_io = None
        img_dct = {}
        wcses = None

        fits_f = None
        try:
            self.logger.info(f'opening file {path}')
            fits_f = FITS(path)

            for idx, hdu in enumerate(fits_f):
                hdu_name = hdu.get_extname()
                if hdu_name == 'PRIMARY':
                    # anything we need to do here?  Something with the header?
                    pass

                elif hdu_name.startswith('CAM'):
                    # load camera image
                    name = hdu_name.strip()
                    cam_num = int(name[-1])

                    #imname = fname + f'[{name}]'
                    imname = f'{name}'
                    image = AstroImage(logger=self.logger,
                                       ioclass=FitsioFileHandler,
                                       wcsclass=AstropyWCS)
                    image.load_hdu(hdu)
                    image.set(name=imname)
                    data = image.get_data()
                    if len(data) == 0:
                        # <-- empty data area--possibly dead camera
                        continue

                    # make a WCS for the image if it doesn't have one
                    if (self._in_gen2 and (is_raw or image.wcs is None or
                                           image.wcs.wcs is None)):
                        # create wcses from Kawanomoto-san's module
                        if wcses is None:
                            wcses = self.make_WCSes()
                        header = image.get_header()
                        wcs = wcses[cam_num - 1]
                        image.update_keywords(wcs.to_header(relax=True))

                    if name == 'CAM1' and set_1k:
                        image.set(path=path)
                        self.fv.gui_do(self.set_1k, image)

                    # perform any desired subtractions
                    image = self.quick_data_reduce(image, name)
                    #image.set(path=f"{path}[{idx}]")

                    # for later updating the images in the grid of viewers
                    img_dct[name] = image

                elif hdu_name == 'detected_objects':
                    self.logger.info('reading do table')
                    self.tbl_do = hdu.read()

                elif hdu_name == 'guide_objects':
                    self.logger.info('reading go table')
                    self.tbl_go = hdu.read()

                elif hdu_name == 'identified_objects':
                    self.logger.info('reading io table')
                    self.tbl_io = hdu.read()

                else:
                    self.logger.info("Unrecognized HDU: name='{}'".format(hdu_name))
            self.logger.info('read all HDUs')

            # determine max and min magnitude
            if self.tbl_go is not None:
                mags = self.tbl_go['mag']
                #self.mag_max = np.max(mags)
                #self.mag_min = np.min(mags)

            self.logger.info('determined minmax')
            end_time = time.time()
            self.logger.info("file processing time %.4f sec" % (end_time - start_time))

            # update the grid of all viewers with new images
            self.fv.gui_do_oneshot('pfsag_update', self.update_grid, img_dct)

        except Exception as e:
            self.logger.error("Failed to process image: {}".format(e),
                              exc_info=True)

        finally:
            fits_f = None

    def make_WCSes(self):
        self.logger.info('fetching status ...')
        ra, dec, pa = self.sc.fetch_list(['STATS.RA_DEG', 'STATS.DEC_DEG',
                                          'FITS.SBR.INST_PA'])
        pa = float(pa)

        self.logger.info('making wcses ...')
        wcses = pfswcs.agcwcs_sip(ra, dec, pa)
        self.logger.info('returning wcses ...')
        return wcses

    def remove_file(self, path):
        try:
            #os.chmod(path, 0o660)
            #os.remove(path)
            pass
        except IOError as e:
            self.logger.error(f"Error removing file {path}: {e}")

    def read_calib(self, path, dct):
        """Read a calibration file and store the slices in dictionary `dct`.
        """
        try:
            with FITS(path) as fits_f:

                for idx, hdu in enumerate(fits_f):
                    hdu_name = hdu.get_extname()
                    if hdu_name.startswith('CAM'):
                        # load camera image
                        dct[hdu_name] = hdu.read()
            self.logger.info("calibration file read successfully")

        except Exception as e:
            self.logger.error("Failed to open calib file '{}': {}".format(path, e),
                              exc_info=True)
        finally:
            fits_f = None

    def clear_stars(self, redraw=True):
        # delete previously plotted objects
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                canvas = viewer.get_canvas()
                canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_go'),
                                             redraw=redraw)
                canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_do'),
                                             redraw=redraw)
                canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_io'),
                                             redraw=redraw)

        # Update PFS_FOV channel
        if self.settings.get('plot_fov', False):
            channel = self.fv.get_channel_on_demand(self.fov_chname)
            viewer = channel.viewer
            fov_canvas = viewer.get_canvas()
            fov_canvas.delete_objects_by_tag(fov_canvas.get_tags_by_tag_pfx('_io'),
                                             redraw=redraw)

    def plot_stars(self):
        try:
            self.clear_stars(redraw=False)
        except Exception as e:
            self.logger.error("Error clearing stars: {e}", exc_info=True)

        plot_fov = self.settings.get('plot_fov', False)
        if plot_fov:
            channel = self.fv.get_channel_on_demand(self.fov_chname)
            viewer = channel.viewer
            fov_canvas = viewer.get_canvas()

        # self.cbar.set_range(self.mag_min, self.mag_max)

        if self.tbl_io is None:
            return
        mags = self.tbl_go['mag']

        do_used = set(self.tbl_io['detected_object_id'])
        go_used = set(self.tbl_io['guide_object_id'])
        #go_not_used = set(range(0, len(self.tbl_go))) - go_used
        go_not_used = set(range(0, len(self.tbl_go)))
        do_not_used = set(range(0, len(self.tbl_do))) - do_used

        radius = 15

        if self.settings.get('plot_guide_stars', False):
            # plot guide objects that are not identified
            color = self.settings.get('guide_color', 'cyan')
            for go_idx in go_not_used:
                go_row = self.tbl_go[go_idx]
                cam_num = go_row['camera_id']
                pos_x, pos_y = (go_row['guide_object_xdet'],
		                go_row['guide_object_ydet'])
                cam_id = 'CAM{}'.format(cam_num + 1)
                if self.fv.has_channel(cam_id):
                    p = self.dc.Point(pos_x, pos_y, radius=radius,
                                      style='diamond',
                                      color=color, linewidth=2)
                    channel = self.fv.get_channel(cam_id)
                    viewer = channel.fitsimage
                    canvas = viewer.get_canvas()
                    canvas.add(p, tag=f'_go{go_idx}', redraw=False)

        if self.settings.get('plot_detected_not_identified', False):
            # plot detected objects that are not identified
            color = self.settings.get('detected_color', 'yellow')
            for do_idx in do_not_used:
                do_row = self.tbl_do[do_idx]
                cam_num = do_row['camera_id']
                # camera indexes are now 0-based, while HDUs are numbered from 1
                cam_id = 'CAM{}'.format(cam_num + 1)
                ctr_x, ctr_y = do_row['centroid_x'], do_row['centroid_y']
                p = self.dc.Point(ctr_x, ctr_y, radius=radius, style='hexagon',
                                  color=color, linewidth=2)

                if self.fv.has_channel(cam_id):
                    channel = self.fv.get_channel(cam_id)
                    viewer = channel.fitsimage
                    canvas = viewer.get_canvas()
                    canvas.add(p, tag=f'_do{do_idx}', redraw=False)

        # plot identified objects
        for io_idx, io_row in enumerate(self.tbl_io):
            do_row = self.tbl_do[io_row['detected_object_id']]
            cam_num = do_row['camera_id']
            # camera indexes are now 0-based, while HDUs are numbered from 1
            cam_id = 'CAM{}'.format(cam_num + 1)
            ctr_x, ctr_y = do_row['centroid_x'], do_row['centroid_y']

            _go_row_num = int(io_row['guide_object_id'])
            mag = mags[_go_row_num]

            # skip stars that fall outside the selected area
            # if not (self.mag_min <= mag <= self.mag_max):
            #     continue

            # add circle for detected position
            #color = self.get_color(mag)
            color = self.settings.get('identified_color', 'orangered')
            objs = []
            fov_objs = []

            if self.settings.get('plot_identified_stars', False):
                c = self.dc.Circle(ctr_x, ctr_y, radius,
                                   color=color, linewidth=2)
                p = self.dc.Point(ctr_x, ctr_y, radius, style='plus',
                                  color=color, linewidth=2)
                objs.extend([c, p])

                image = self.img_dct[cam_id]

                if plot_fov:
                    fov_ctr_x, fov_ctr_y = self.get_fov_xy(cam_num, ctr_x, ctr_y)
                    c = self.dc.Circle(fov_ctr_x, fov_ctr_y, radius,
                                       color=color, linewidth=2)
                    p = self.dc.Point(fov_ctr_x, fov_ctr_y, radius, style='plus',
                                      color=color, linewidth=2)
                    fov_objs.extend([c, p])

                if self.settings.get('plot_offsets', False):

                    gde_x, gde_y = (io_row['guide_object_xdet'],
                                    io_row['guide_object_ydet'])
                    # go_row = self.tbl_go[_go_row_num]
                    # gde_x, gde_y = (go_row['guide_object_xdet'],
		    #                 go_row['guide_object_ydet'])

                    error = np.sqrt((gde_y - ctr_y) ** 2 + (gde_x - ctr_x) ** 2)
                    # scale the error for better visibility
                    err_long = error * self.error_scale

                    theta_rad = np.arctan2(gde_y - ctr_y, gde_x - ctr_x)
                    long_x, long_y = (ctr_x + err_long * np.cos(theta_rad),
                                      ctr_y + err_long * np.sin(theta_rad))

                    c = self.dc.Circle(gde_x, gde_y, radius,
                                       color=color, linestyle='dash',
                                       linewidth=2)
                    #l = self.dc.Line(ctr_x, ctr_y, gde_x, gde_y,
                    l = self.dc.Line(ctr_x, ctr_y, long_x, long_y,
                                     color=color, linestyle='solid',
                                     linewidth=2, arrow='end')
                    objs.extend([c, l])

                    if plot_fov:
                        fov_gde_x, fov_gde_y = self.get_fov_xy(cam_num, gde_x, gde_y)

                        theta_rad = np.arctan2(fov_gde_y - fov_ctr_y, fov_gde_x - fov_ctr_x)
                        fov_long_x, fov_long_y = (fov_ctr_x + err_long * np.cos(theta_rad),
                                                  fov_ctr_y + err_long * np.sin(theta_rad))

                        c = self.dc.Circle(fov_gde_x, fov_gde_y, radius,
                                           color=color, linestyle='dash',
                                           linewidth=2)
                        #l = self.dc.Line(fov_ctr_x, fov_ctr_y, fov_gde_x, fov_gde_y,
                        l = self.dc.Line(fov_ctr_x, fov_ctr_y, fov_long_x, fov_long_y,
                                         color=color, linestyle='solid',
                                         linewidth=2, arrow='end')
                        fov_objs.extend([c, l])

            if len(objs) > 0:
                if self.fv.has_channel(cam_id):
                    channel = self.fv.get_channel(cam_id)
                    viewer = channel.viewer
                    canvas = viewer.get_canvas()

                    canvas.add(self.dc.CompoundObject(*objs),
                               tag=f'_io{io_idx}', redraw=False)

            if len(fov_objs) > 0 and plot_fov:
                if self.fv.has_channel(self.fov_chname):
                    fov_canvas.add(self.dc.CompoundObject(*fov_objs),
                                   tag=f'_io{io_idx}', redraw=False)

        # update all canvases
        if self.fv.has_channel(self.fov_chname) and plot_fov:
            fov_canvas.update_canvas(whence=3)

        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                canvas = viewer.get_canvas()
                canvas.update_canvas(whence=3)

    # def get_color(self, mag):
    #     # calculate range of values
    #     rng = float(self.mag_max - self.mag_min)

    #     # clip magnitude to the range we have defined
    #     mag = np.clip(mag, self.mag_min, self.mag_max)

    #     if rng != 0.0:
    #         point = float(mag - self.mag_min) / rng
    #     else:
    #         point = 1.0

    #     # sanity check: clip to 8-bit color range
    #     point = int(np.clip(point * 255.0, 0, 255))
    #     #point = int(point * 255.0)

    #     # Apply colormap.
    #     rgbmap = self.cbar.get_rgbmap()
    #     (r, g, b) = rgbmap.get_rgbval(point)
    #     r = float(r) / 255.0
    #     g = float(g) / 255.0
    #     b = float(b) / 255.0
    #     return (r, g, b)

    def replot_stars(self, rgbmap):
        # not necessary if link=True when creating ColorBar object
        #self.cbar.cbar_view.redraw()
        self.plot_stars()

    def toggle_plot_fov_cb(self, w, tf):
        self.settings.set(plot_fov=tf)

    def toggle_plot_identified_stars_cb(self, w, tf):
        self.settings.set(plot_identified_stars=tf)
        self.plot_stars()

    def toggle_plot_offsets_cb(self, w, tf):
        self.settings.set(plot_offsets=tf)
        self.plot_stars()

    def toggle_plot_guide_stars_cb(self, w, tf):
        self.settings.set(plot_guide_stars=tf)
        self.plot_stars()

    def toggle_plot_ni_detected_stars_cb(self, w, tf):
        self.settings.set(plot_detected_not_identified=tf)
        self.plot_stars()

    def subtract_bias_cb(self, w, tf):
        self.settings.set(subtract_bias=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def subtract_bg_cb(self, w, tf):
        self.settings.set(subtract_background=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    # def subtract_dark_cb(self, w, tf):
    #     self.settings.set(subtract_dark=tf)
    #     if self.current_file is not None:
    #         self.fv.nongui_do(self.process_file, self.current_file)

    # def divide_flat_cb(self, w, tf):
    #     self.settings.set(divide_flat=tf)
    #     if self.current_file is not None:
    #         self.fv.nongui_do(self.process_file, self.current_file)

    def set_flats_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.flat)

    def set_darks_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.dark)

    # def set_minmax_cb(self, *args):
    #     minval = float(self.w.minf.get_text().strip())
    #     self.mag_min = minval
    #     maxval = float(self.w.maxf.get_text().strip())
    #     self.mag_max = maxval

    #     self.plot_stars()

    def save_file_cb(self, w):
        _dir, filename = os.path.split(self.current_file)
        dst_path = os.path.join(self.save_dir, filename)
        res = os.system("cp {} {}".format(self.current_file, dst_path))
        if res == 0 and os.path.exists(dst_path):
            channel = self.fv.get_channel_on_demand('PFS_SAVE')
            image = AstroImage(logger=self.logger)
            image.load_file(dst_path)
            channel.add_image(image)

    def redo(self, channel, image):
        """This is called when a new image arrives or the data in the
        existing image changes.
        """
        ## if image is not None and image.io.fileinfo is not None:
        ##     self.process_image(image)
        pass

    def watch_loop(self, ev_quit):

        data_dir = self.settings.get('data_directory', '.')
        i = inotify.adapters.Inotify()
        i.add_watch(data_dir)

        while not ev_quit.is_set():
            events = list(i.event_gen(yield_nones=False, timeout_s=1.0))
            loop_tot = []
            fits_tot = []

            for event in events:
                if event is not None:
                    if ('IN_MOVED_TO' in event[1] or
                        'IN_CLOSE_WRITE' in event[1]):
                        (header, type_names, watch_path, filename) = event
                        filepath = os.path.join(watch_path, filename)
                        filedir, filename = os.path.split(filepath)
                        loop_tot.append(filename)
                        # sanity check--this is a FITS file, right?
                        if not filename.endswith('.fits'):
                            continue

                        # save filepath for last received one of each type
                        if filename.startswith('_'):
                            self.last.setvals(raw=filepath,
                                              time_raw=time.time())
                            if self.mode == 'processed':
                                continue
                        else:
                            self.last.setvals(processed=filepath,
                                              time_processed=time.time())
                            if self.mode == 'raw':
                                continue

                        start_time = time.time()
                        if start_time < self.last.image_time + self.rate_limit:
                            self.logger.info(f"skipping file '{filename}' for rate limit")
                            #self.remove_file(filepath)
                            continue

                        self.last.image_time = start_time
                        self.logger.debug(f"new file detected: '{filename}'")
                        if self.pause_flag:
                            self.logger.debug("plugin is paused, skipping new file")
                            continue

                        fits_tot.append(filename)

                        self.fv.nongui_do(self.process_file, filepath,
                                          set_1k=True)

            self.logger.debug("---------")
            self.logger.debug("{} files: {}".format(len(loop_tot), loop_tot))
            self.logger.debug("{} fits: {}".format(len(fits_tot), fits_tot))

        i.remove_watch(data_dir)

    # def set_cmap_cb(self, w, idx):
    #     # Get colormap
    #     name = w.get_text()
    #     cm = cmap.get_cmap(name)
    #     self.cbar.set_cmap(cm)
    #     self.plot_stars()

    # def set_imap_cb(self, w, idx):
    #     # Get intensity map
    #     name = w.get_text()
    #     im = imap.get_imap(name)
    #     self.cbar.set_imap(im)
    #     self.plot_stars()

    # def set_field_cb(self, w, idx):
    #     # Get field name
    #     self.field = w.get_text().lower()
    #     self.plot_stars()

    def set_error_scale_cb(self, w, val):
        self.error_scale = 1.0 + val * 0.25
        self.w.error_scale.set_text("{:.2f}".format(self.error_scale))
        self.plot_stars()

    def set_error_scale2_cb(self, w):
        self.error_scale = float(w.get_text())
        i = int(max(0.0, min(50.0, (self.error_scale - 1) / 0.25)))
        self.w.error_scaling.set_value(i)
        self.plot_stars()

    def pause_cb(self, w, tf):
        self.pause_flag = tf

    def auto_orient(self):
        auto_orient = self.settings.get('auto_orient', False)
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                with viewer.suppress_redraw:
                    viewer.transform(False, False, False)
                    viewer.rotate(0.0)
                    if auto_orient:
                        rot_ang_deg = self.rot_angles[cam_num - 1]
                        viewer.rotate(rot_ang_deg)
                        viewer.transform(False, True, False)

    def auto_orient_cb(self, w, tf):
        self.settings.set(auto_orient=tf)
        self.auto_orient()

    def set_rate_limit_cb(self, w, val):
        self.rate_limit = val

    def set_mode_cb(self, w, tf, mode):
        if not tf:
            return
        old_mode = self.mode
        self.mode = mode

        if old_mode != mode:
            # process the last known received file of the kind we are
            # switching to
            filepath = self.last.raw if mode == 'raw' else self.last.processed
            if filepath is None or not os.path.exists(filepath):
                return
            self.fv.nongui_do(self.process_file, filepath, set_1k=True)

    def pan_cam_cb(self, w, cam_num):
        channel = self.fv.get_channel(self.fov_chname)
        viewer = channel.viewer
        if cam_num not in self._fov_coords:
            viewer.onscreen_message(f"CAM{cam_id+1} image not available",
                                    delay=1.0)
            return

        base_x, base_y, sub_x, sub_y = self._fov_coords[cam_num]
        with viewer.suppress_redraw:
            viewer.set_pan(base_x, base_y, coord='data')
            viewer.scale_to(0.65, 0.65)

    def pan_zoom_fit_cb(self, w):
        channel = self.fv.get_channel(self.fov_chname)
        viewer = channel.viewer
        viewer.zoom_fit()

    def connect_status(self):
        self.sc = StatusClient(host=self.config_d['status_client_host'],
                               username=self.config_d['status_client_user'],
                               password=self.config_d['status_client_pass'])
        self.sc.connect()

    def __str__(self):
        return 'pfs_ag'
