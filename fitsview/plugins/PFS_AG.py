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
from ginga.util.wcsmod.wcs_astropy import AstropyWCS
# pip install inotify
import inotify.adapters
import yaml

# ginga
from ginga import trcalc, cmap, imap
from ginga.gw import ColorBar, Widgets, Viewers
from ginga.AstroImage import AstroImage
from ginga.util.io_fits import FitsioFileHandler
from ginga.util.wcsmod.wcs_astropy import AstropyWCS
from ginga.util import wcs, dp
from ginga.util.mosaic import CanvasMosaicer
from ginga import GingaPlugin

# g2cam
from g2cam.status.client import StatusClient

# local
from fitsview.util import pfswcs

is_summit = True


class PFS_AG(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(PFS_AG, self).__init__(fv)

        # get PFSAG preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_PFS_AG')
        self.settings.add_defaults(plot_fov=True,
                                   collage_method='simple',
                                   plot_guide_stars=False,
                                   plot_offsets=False,
                                   subtract_bias=False,
                                   subtract_background=False,
                                   color_map='rainbow3',
                                   intensity_map='ramp',
                                   channel_name='PFS_1K',
                                   rate_limit=5.0,
                                   data_directory='.',
                                   save_directory=tempfile.gettempdir(),
                                   auto_orient=False)
        self.settings.load(onError='silent')

        self.chname = self.settings.get('channel_name')
        self.wsname = 'PFS_AG_CAMS'
        self.wstype = 'grid'
        self.inspace = 'channels' if not is_summit else 'sub1' # 'top_level'
        self.fov_chname = 'PFS_FOV'
        self.fov_inspace = 'channels' if not is_summit else 'sub2'

        self._wd = 300
        self._ht = 300
        self.mag_max = 20.0
        self.mag_min = 12.0
        self.current_file = None
        self.ev_quit = threading.Event()
        self.dark = dict()
        self.flat = dict()
        self.cmap_names = list(cmap.get_names())
        self.cmap = cmap.get_cmap(self.settings.get('color_map'))
        self.imap_names = list(imap.get_names())
        self.imap = imap.get_imap(self.settings.get('intensity_map'))
        self.field_names = ['mag']
        self.field = 'mag'
        self.last_image_time = time.time()
        self.pause_flag = False
        self.rate_limit = self.settings.get('rate_limit', 5.0)
        self.error_scale = 10
        self.save_dir = self.settings.get('save_directory', '/tmp')
        self.mode = 'processed'
        self.guide_count = 0

        # hold tables of detected objs, guide objs and identified objs
        self.tbl_do = None
        self.tbl_go = None
        self.tbl_io = None
        self.img_dct = {}

        self.viewer = dict()
        self.dc = fv.get_draw_classes()

        self.mosaicer = CanvasMosaicer(self.logger)
        t_ = self.mosaicer.get_settings()
        collage_method = self.settings['collage_method']
        t_.set(annotate_images=True, match_bg=False,
               center_image=False, collage_method=collage_method)

        if not is_summit:
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
                lbl = self.dc.Text(20, 30, chname, color='skyblue',
                                   fontsize=12, font='sans',
                                   coord='window')
                lbl.crdmap = viewer.get_coordmap('window')
                canvas = viewer.get_canvas()
                canvas.add(lbl, redraw=False)
                settings = viewer.get_settings()
                settings.set(enter_focus=True) #, channel_follows_focus=True)
                #viewer.add_callback('focus', self.focus_cb)

        # create "big picture" FOV channel
        chname = self.fov_chname
        settings = prefs.create_category(f'channel_{chname}')
        settings.set(numImages=1, genthumb=False, raisenew=False)
        settings.load(onError='silent')
        channel = self.fv.add_channel(chname, settings=settings,
                                      workspace=self.fov_inspace)
        viewer = channel.viewer
        #canvas = viewer.get_canvas()
        canvas = self.dc.DrawingCanvas()
        ## canvas.add_draw_mode('zoom', key=self.zoom_fov_cb)
        ## canvas.register_for_cursor_drawing(viewer)
        ## viewer.get_canvas().add(canvas, redraw=False)
        ## canvas.set_draw_mode('zoom')
        canvas.set_surface(viewer)
        canvas.ui_set_active(True)
        self.fov_canvas = canvas

        captions = [('Plot FOV', 'checkbutton',
                     "Method:", 'label', 'method', 'combobox'),
                    ('Plot guide stars', 'checkbutton',
                     'Plot offsets', 'checkbutton',
                     'Auto orient', 'checkbutton'),
                    ('Subtract bias', 'checkbutton'),
                    ('Subtract dark', 'checkbutton', 'Darks:', 'llabel',
                     'dark_frames', 'entryset'),
                    ('Subtract background', 'checkbutton'),
                    ('Divide flat', 'checkbutton', 'Flats:', 'llabel',
                     'flat_frames', 'entryset'),
                    ('Pause', 'togglebutton', 'Rate Limit', 'spinfloat',
                     'Save', 'button'),
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        tf = self.settings.get('plot_fov', False)
        b.plot_fov.set_state(tf)
        b.plot_fov.set_tooltip("Plot images/stars in PFS FOV")
        b.plot_fov.add_callback('activated', self.toggle_plot_fov_cb)

        combobox = b.method
        options = ['simple', 'warp']
        for name in options:
            combobox.append_text(name)
        method = self.settings.get('collage_method', 'simple')
        combobox.set_text(method)
        combobox.add_callback('activated', self.set_collage_method_cb)
        combobox.set_tooltip("Choose collage method: %s" % ','.join(options))

        tf = self.settings.get('plot_guide_stars', False)
        b.plot_guide_stars.set_state(tf)
        b.plot_guide_stars.set_tooltip("Plot the guide stars")
        b.plot_guide_stars.add_callback('activated', self.toggle_plot_guide_stars_cb)

        tf = self.settings.get('plot_offsets', False)
        b.plot_offsets.set_state(tf)
        b.plot_offsets.set_tooltip("Plot the offsets from the guide stars")
        b.plot_offsets.add_callback('activated', self.toggle_plot_offsets_cb)

        tf = self.settings.get('auto_orient', False)
        b.auto_orient.set_state(tf)
        b.auto_orient.set_tooltip("Auto rotate images to orient by N")
        b.auto_orient.add_callback('activated', self.auto_orient_cb)

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

        tf = self.settings.get('subtract_dark', False)
        b.subtract_dark.set_state(tf)
        b.subtract_dark.add_callback('activated', self.subtract_dark_cb)
        b.subtract_dark.set_tooltip("Subtract dark from frame")
        b.dark_frames.add_callback('activated', self.set_darks_cb)
        b.dark_frames.set_tooltip("Enter a file that contains dark frames")
        #b.subtract_dark.set_enabled(False)
        top.add_widget(w, stretch=0)

        tf = self.settings.get('divide_flat', False)
        b.divide_flat.set_state(tf)
        b.divide_flat.add_callback('activated', self.divide_flat_cb)
        b.divide_flat.set_tooltip("Divide frame by flat")
        #b.divide_flat.set_enabled(False)
        b.flat_frames.add_callback('activated', self.set_flats_cb)
        b.flat_frames.set_tooltip("Enter a file that contains flat frames")

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

        # add CAM pan buttons
        hbox = Widgets.HBox()
        for i in range(1, 7):
            w = Widgets.Button(f"CAM{i}")
            hbox.add_widget(w, stretch=1)
            w.add_callback('activated', self.pan_cam_cb, i)
            w.set_tooltip(f"Pan to CAM{i} in the PFS_FOV viewer")
        top.add_widget(hbox, stretch=0)

        # stretch spacer
        top.add_widget(Widgets.Label(''), stretch=1)

        captions = [('error_scaling', 'hslider'),
                    ('Set minmax', 'button', 'minf', 'entry', 'maxf', 'entry'),
                    ('color_map', 'combobox', 'intensity_map', 'combobox',
                     'field', 'combobox')
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.error_scaling.set_limits(1, 100, incr_value=5)
        b.error_scaling.set_value(self.error_scale)
        b.error_scaling.add_callback('value-changed', self.set_error_scale_cb)
        b.error_scaling.set_tooltip("Change the error scaling")

        b.minf.set_text(str(self.mag_min))
        b.minf.add_callback('activated', self.set_minmax_cb)
        b.maxf.set_text(str(self.mag_max))
        b.maxf.add_callback('activated', self.set_minmax_cb)
        b.set_minmax.add_callback('activated', self.set_minmax_cb)
        b.set_minmax.set_tooltip("Set min/max values for color bar")

        combobox = b.color_map
        for name in self.cmap_names:
            combobox.append_text(name)
        cmap_name = 'rainbow3'
        try:
            index = self.cmap_names.index(cmap_name)
        except Exception:
            index = self.cmap_names.index('gray')
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_cmap_cb)

        combobox = b.intensity_map
        for name in self.imap_names:
            combobox.append_text(name)
        imap_name = 'ramp'
        try:
            index = self.imap_names.index(imap_name)
        except Exception:
            index = self.imap_names.index('ramp')
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_imap_cb)

        combobox = b.field
        for name in self.field_names:
            combobox.append_text(name)
        field_name = 'mag'
        try:
            index = self.field_names.index(field_name)
        except Exception:
            index = self.field_names.index('mag')
        combobox.set_index(index)
        combobox.add_callback('activated', self.set_field_cb)

        top.add_widget(w, stretch=0)

        self.cbar = ColorBar.ColorBar(self.logger, link=True)
        self.cbar.set_cmap(self.cmap)
        self.cbar.set_imap(self.imap)
        rgbmap = self.cbar.get_rgbmap()
        rgbmap.add_callback('changed', self.replot_stars)
        # hack to set font size of this color bar
        self.cbar.cbar.fontsize = 8

        cbar_w = self.cbar.get_widget()
        cbar_w.resize(-1, 32)
        top.add_widget(cbar_w, stretch=0)

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
        if is_summit:
            self.fv.nongui_do(self.watch_loop, self.ev_quit)

        self.connect_status()

    def stop(self):
        self.ev_quit.set()
        self.gui_up = False

    def focus_cb(self, viewer, tf):
        if not tf:
            return
        chname = self.fv.get_channel_name(viewer)
        self.fv.change_channel(chname, raisew=True)

    def incoming_data_cb(self, fv, chname, image, info):
        if chname != 'Image':
            return

        #header = image.get_header()

        # add image to obslog
        #self.fv.gui_do(self.add_to_obslog, header, image)

        self.fv.nongui_do(self.process_image, image)

    def quick_data_reduce(self, image, name):
        data = image.get_data()
        wd, ht = image.get_size()[:2]
        ovsc_x, ovsc_y = 24, 9

        if self.settings.get('subtract_bias', False):
            # overscans - 24 pixels at the beginning and end of each row (x)
            # and first 9 rows (y).
            overscan = np.concatenate([np.ravel(data[0:ovsc_x, 0:ht]),
                                       np.ravel(data[wd - ovsc_x:wd, 0:ht]),
                                       np.ravel(data[ovsc_x:wd - ovsc_x, 0:ovsc_y])])
            med = np.nanmedian(overscan)
            data = data - med

        if self.settings.get('subtract_background', False):
            med = np.nanmedian(data[ovsc_x:wd - ovsc_x, ovsc_y:ht])
            data = data - med

        if self.settings.get('subtract_dark', False):
            if len(self.dark) == 0:
                self.logger.error("No darks are loaded")
            else:
                data = data - self.dark[name]

        if self.settings.get('divide_flat', False):
            if len(self.dark) == 0:
                self.logger.error("No flats are loaded")
            else:
                data = data / self.flat[name]

        new_img = AstroImage(data_np=data)
        new_img.update_keywords(image.get_header())
        new_img.set(name=image.get('name'), nothumb=True)
        #new_img.set(name='PFS' + str(time.time()), nothumb=True)
        return new_img

    def process_image(self, image):
        path = image.get('path', None)
        if path is None:
            return

        self.process_file(path, set_1k=True)

    def orient(self, viewer):
        if self.settings.get('auto_orient', False):
            bd = viewer.get_bindings()
            mode = bd.get_mode_obj('rotate')
            mode._orient(viewer, righthand=False, msg=False)

    def update_grid(self, img_dct):
        self.img_dct = img_dct
        # NOTE: assumes images come in the order CAM1 .. CAM6
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
                    #canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_io'),
                    #                             redraw=False)
                    #viewer.set_image(image)
                    channel.add_image(image)
                    self.orient(viewer)
                else:
                    viewer.clear()

        self.fv.update_pending()

        if is_summit:
            # increment the guide count
            self.guide_count += 1
            stat_d = {'VGW.PFS.AG.COUNT': self.guide_count}
            self.fv.nongui_do(self.fv.call_global_plugin_method,
                              'Gen2Int', 'store', [stat_d], {})

        if self.settings.get('plot_fov', False):
            images = list(img_dct.values())
            #ref_image = self.get_center(images)
            #ref_image.set(tag='CAM0')
            ref_image = images[0]
            self.ref_image = ref_image

            channel = self.fv.get_channel_on_demand(self.fov_chname)
            viewer = channel.viewer
            canvas = viewer.get_canvas()
            ## canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('CAM'),
            ##                              redraw=False)
            ## canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_io'),
            ##                              redraw=False)

            self.mosaicer.reset()
            with viewer.suppress_redraw:
                self.mosaicer.mosaic(viewer, images, canvas=canvas)
                self.orient(viewer)
                viewer.redraw(whence=0)

        if self.tbl_go is not None:
            self.plot_stars()

    def set_1k(self, image):
        self.fv.assert_gui_thread()

        channel = self.fv.get_channel_on_demand(self.chname)
        channel.add_image(image)

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

                    self.logger.info(f'{fname}, {is_raw}')
                    # make a WCS for the image if it doesn't have one
                    if (is_summit and (is_raw or image.wcs is None or
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

                    # update the image in the channel viewer
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

    def get_center(self, images):
        coords = []
        for im in images:
            wd, ht = im.get_size()
            coords.append((im.pixtoradec(wd * 0.5, ht * 0.5)))
        coords = np.array(coords)
        ra_deg, dec_deg = np.mean(coords, axis=0)

        fov_deg = 0.001
        px_scale = 4.0694108718577e-05    # measured
        rot_deg = 0.0

        img_ctr = dp.create_blank_image(ra_deg, dec_deg, fov_deg,
                                        px_scale, rot_deg, dtype=np.uint)
        return img_ctr


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
            with fits.open(path, 'readonly') as fits_f:

                for idx, hdu in enumerate(fits_f):
                    hdu_name = hdu.get_extname()
                    if hdu_name.startswith('CAM'):
                        # load camera image
                        dct[hdu_name] = hdu.read()

        except Exception as e:
            self.logger.error("Failed to open calib file '{}': {}".format(path, e),
                              exc_info=True)
        finally:
            fits_f = None

    def plot_stars(self):
        self.cbar.set_range(self.mag_min, self.mag_max)

        if self.tbl_io is None:
            return
        mags = self.tbl_go['mag']

        # delete previously plotted objects
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                canvas = viewer.get_canvas()
                canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_io'),
                                             redraw=False)

        # Update PFS_FOV channel
        if self.settings.get('plot_fov', False):
            channel = self.fv.get_channel_on_demand(self.fov_chname)
            viewer = channel.viewer
            fov_canvas = viewer.get_canvas()
            fov_canvas.delete_objects_by_tag(fov_canvas.get_tags_by_tag_pfx('_io'),
                                             redraw=False)
        ref_image = self.ref_image

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
            if not (self.mag_min <= mag <= self.mag_max):
                continue

            # add circle for detected position
            color = self.get_color(mag)
            radius = 15
            objs = []
            fov_objs = []

            if self.settings.get('plot_guide_stars', False):
                c = self.dc.Circle(ctr_x, ctr_y, radius,
                                   color=color, linewidth=2)
                p = self.dc.Point(ctr_x, ctr_y, radius, style='plus',
                                  color=color, linewidth=2)
                objs.extend([c, p])

                image = self.img_dct[cam_id]

                if self.settings.get('plot_fov', False):
                    ra_ctr, dec_ctr = image.pixtoradec(ctr_x, ctr_y)
                    fov_ctr_x, fov_ctr_y = ref_image.radectopix(ra_ctr, dec_ctr)

                    c = self.dc.Circle(fov_ctr_x, fov_ctr_y, radius,
                                       color=color, linewidth=2)
                    p = self.dc.Point(fov_ctr_x, fov_ctr_y, radius, style='plus',
                                      color=color, linewidth=2)
                    fov_objs.extend([c, p])

                if self.settings.get('plot_offsets', False):

                    gde_x, gde_y = (io_row['guide_object_xdet'],
                                    io_row['guide_object_ydet'])

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

                    if self.settings.get('plot_fov', False):
                        ra_gde, dec_gde = image.pixtoradec(gde_x, gde_y)
                        fov_gde_x, fov_gde_y = ref_image.radectopix(ra_gde, dec_gde)
                        ## ra_deg = self.tbl_go['ra'][_go_row_num]
                        ## dec_deg = self.tbl_go['dec'][_go_row_num]
                        ## fov_gde_x, fov_gde_y = ref_image.radectopix(ra_deg, dec_deg)

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

            if len(fov_objs) > 0 and self.settings.get('plot_fov', False):
                if self.fv.has_channel(self.fov_chname):
                    fov_canvas.add(self.dc.CompoundObject(*fov_objs),
                                   tag=f'_io{io_idx}', redraw=False)

        # update all canvases
        if self.fv.has_channel(self.fov_chname):
            fov_canvas.update_canvas(whence=3)

        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                canvas = viewer.get_canvas()
                canvas.update_canvas(whence=3)

    def get_color(self, mag):

        # calculate range of values
        rng = float(self.mag_max - self.mag_min)

        # clip magnitude to the range we have defined
        mag = np.clip(mag, self.mag_min, self.mag_max)

        if rng != 0.0:
            point = float(mag - self.mag_min) / rng
        else:
            point = 1.0

        # sanity check: clip to 8-bit color range
        point = int(np.clip(point * 255.0, 0, 255))
        #point = int(point * 255.0)

        # Apply colormap.
        rgbmap = self.cbar.get_rgbmap()
        (r, g, b) = rgbmap.get_rgbval(point)
        r = float(r) / 255.0
        g = float(g) / 255.0
        b = float(b) / 255.0
        return (r, g, b)

    def replot_stars(self, rgbmap):
        # not necessary if link=True when creating ColorBar object
        #self.cbar.cbar_view.redraw()
        self.plot_stars()

    def toggle_plot_fov_cb(self, w, tf):
        self.settings.set(plot_fov=tf)

    def toggle_plot_guide_stars_cb(self, w, tf):
        self.settings.set(plot_guide_stars=tf)
        self.plot_stars()

    def toggle_plot_offsets_cb(self, w, tf):
        self.settings.set(plot_offsets=tf)
        self.plot_stars()

    def subtract_bias_cb(self, w, tf):
        self.settings.set(subtract_bias=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def subtract_bg_cb(self, w, tf):
        self.settings.set(subtract_background=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def subtract_dark_cb(self, w, tf):
        self.settings.set(subtract_dark=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def divide_flat_cb(self, w, tf):
        self.settings.set(divide_flat=tf)
        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def set_flats_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.flat)

    def set_darks_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.dark)

    def set_minmax_cb(self, *args):
        minval = float(self.w.minf.get_text().strip())
        self.mag_min = minval
        maxval = float(self.w.maxf.get_text().strip())
        self.mag_max = maxval

        self.plot_stars()

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

                        # TEMP: choose which type to accept
                        if self.mode == 'processed' and filename.startswith('_'):
                            continue
                        if self.mode == 'raw' and not filename.startswith('_'):
                            continue

                        start_time = time.time()
                        if start_time < self.last_image_time + self.rate_limit:
                            self.logger.info(f"skipping file '{filename}' for rate limit")
                            self.remove_file(filepath)
                            continue
                        self.last_image_time = start_time

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

    def set_cmap_cb(self, w, idx):
        # Get colormap
        name = w.get_text()
        cm = cmap.get_cmap(name)
        self.cbar.set_cmap(cm)
        self.plot_stars()

    def set_imap_cb(self, w, idx):
        # Get intensity map
        name = w.get_text()
        im = imap.get_imap(name)
        self.cbar.set_imap(im)
        self.plot_stars()

    def set_field_cb(self, w, idx):
        # Get field name
        self.field = w.get_text().lower()
        self.plot_stars()

    def set_error_scale_cb(self, w, val):
        self.error_scale = val
        self.plot_stars()

    def zoom_fov_cb(self, canvas, event, *args):
        if event.key in ['1', '2', '3', '4', '5', '6']:
            cam = 'CAM%d' % event.key
            return True
        return False

    def pause_cb(self, w, tf):
        self.pause_flag = tf

    def auto_orient_cb(self, w, tf):
        self.settings.set(auto_orient=tf)
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            if self.fv.has_channel(cam_id):
                channel = self.fv.get_channel(cam_id)
                viewer = channel.viewer
                if tf:
                    self.orient(viewer)
                else:
                    viewer.transform(False, False, False)
                    viewer.rotate(0.0)

        if self.fv.has_channel(self.fov_chname):
            channel = self.fv.get_channel(self.fov_chname)
            viewer = channel.viewer
            if tf:
                self.orient(viewer)
            else:
                viewer.transform(False, False, False)
                viewer.rotate(0.0)

    def set_rate_limit_cb(self, w, val):
        self.rate_limit = val

    def set_mode_cb(self, w, tf, mode):
        if tf:
            self.mode = mode

    def pan_cam_cb(self, w, cam_num):
        cam_id = 'CAM{}'.format(cam_num)
        channel = self.fv.get_channel(self.fov_chname)
        viewer = channel.viewer
        if cam_id not in self.img_dct:
            viewer.onscreen_message(f"{cam_id} image not available",
                                    delay=1.0)
            return

        for name, tag, ra, dec in self.mosaicer.image_list:
            if name == cam_id:
                with viewer.suppress_redraw:
                    viewer.set_pan(ra, dec, coord='wcs')
                    #viewer.zoom_to(-2.0)
                return

    def set_collage_method_cb(self, w, idx):
        method = w.get_text()
        settings = self.mosaicer.get_settings()
        settings.set(collage_method=method)

        # "warp" currently takes longer, so automatically set rate limit
        # to accomodate that
        if method == 'warp':
            self.rate_limit = max(15.0, self.rate_limit)
            self.w.rate_limit.set_value(self.rate_limit)

        if self.current_file is not None:
            self.fv.nongui_do(self.process_file, self.current_file)

    def connect_status(self):
        self.sc = StatusClient(host=self.config_d['status_client_host'],
                               username=self.config_d['status_client_user'],
                               password=self.config_d['status_client_pass'])
        self.sc.connect()

    def __str__(self):
        return 'PFS_AG'
