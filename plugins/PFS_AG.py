#
# PFSAG.py -- PFS AG routines
#
# Eric Jeschke (eric@naoj.org)
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

# 3rd party
import numpy as np
from astropy.io import fits
# pip install inotify
import inotify.adapters

# ginga
from ginga import trcalc, cmap, imap
from ginga.util.io_fits import PyFitsFileHandler
from ginga.gw import ColorBar, Widgets, Viewers
from ginga.AstroImage import AstroImage
from ginga import GingaPlugin


class PFS_AG(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(PFS_AG, self).__init__(fv)

        # get PFSAG preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_PFS_AG')
        self.settings.add_defaults(plot_guide_stars=False,
                                   plot_offsets=False,
                                   subtract_bias=False,
                                   subtract_background=False,
                                   color_map='rainbow3',
                                   intensity_map='ramp',
                                   channel_name='PFS_1K',
                                   rate_limit=5.0,
                                   data_directory='.')
        self.settings.load(onError='silent')

        self.chname = self.settings.get('channel_name')
        self.wsname = 'PFS_AG_CAMS'
        self.wstype = 'grid'
        #self.inspace = 'top level'
        #self.inspace = 'channels'
        self.inspace = 'sub1'

        self._wd = 300
        self._ht = 300
        self.mag_max = 0.0
        self.mag_min = 0.0
        self.current_file = None
        self.ev_quit = threading.Event()
        self.dark = dict()
        self.flat = dict()
        self.cmap_names = list(cmap.get_names())
        self.cmap = cmap.get_cmap(self.settings.get('color_map'))
        self.imap_names = list(imap.get_names())
        self.imap = imap.get_imap(self.settings.get('intensity_map'))
        self.last_image_time = time.time()
        self.pause_flag = False
        self.rate_limit = self.settings.get('rate_limit', 5.0)

        # hold tables of detected objs, guide objs and identified objs
        self.tbl_do = None
        self.tbl_go = None
        self.tbl_io = None

        self.viewer = dict()
        self.processed = set([])
        self.opener = PyFitsFileHandler(self.logger)
        self.dc = fv.get_draw_classes()

        #self.fv.add_callback('add-image', self.incoming_data_cb)
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
                settings = viewer.get_settings()
                settings.set(enter_focus=True) #, channel_follows_focus=True)
                #viewer.add_callback('focus', self.focus_cb)

        captions = [('Plot guide stars', 'checkbutton'),
                    ('Plot offsets', 'checkbutton'),
                    ('Subtract bias', 'checkbutton'),
                    ('Subtract dark', 'checkbutton', 'Darks:', 'llabel',
                     'dark_frames', 'entryset'),
                    ('Subtract background', 'checkbutton'),
                    ('Divide flat', 'checkbutton', 'Flats:', 'llabel',
                     'flat_frames', 'entryset'),
                    ('Pause', 'togglebutton', 'Rate Limit', 'spinfloat'),
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)
        tf = self.settings.get('plot_guide_stars', False)
        b.plot_guide_stars.set_state(tf)
        b.plot_guide_stars.set_tooltip("Plot the guide stars")
        b.plot_guide_stars.add_callback('activated', self.toggle_plot_guide_stars_cb)

        tf = self.settings.get('plot_offsets', False)
        b.plot_offsets.set_state(tf)
        b.plot_offsets.set_tooltip("Plot the offsets from the guide stars")
        b.plot_offsets.add_callback('activated', self.toggle_plot_offsets_cb)

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
        b.rate_limit.set_limits(5.0, 60.0, incr_value=1.0)
        b.rate_limit.set_value(self.rate_limit)
        b.rate_limit.add_callback('value-changed', self.set_rate_limit_cb)
        b.rate_limit.set_tooltip("Set the rate limit for PFS AG files")

        top.add_widget(w, stretch=0)

        # stretch spacer
        top.add_widget(Widgets.Label(''), stretch=1)

        captions = [('color_map', 'combobox',
                     'intensity_map', 'combobox')
                    ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

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
        self.fv.nongui_do(self.watch_loop, self.ev_quit)

    def stop(self):
        self.ev_quit.set()
        self.gui_up = False

    def focus_cb(self, viewer, tf):
        if not tf:
            return
        chname = self.fv.get_channel_name(viewer)
        self.fv.change_channel(chname, raisew=True)

    def incoming_data_cb(self, fv, chname, image, info):
        if chname != self.chname:
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
        #new_img.set(name=image.get('name'))
        new_img.set(name='PFS' + str(time.time()), nothumb=True)
        return new_img

    def process_image(self, image):
        imname = image.get('name', None)
        if imname is None:
            return

        path = image.io.fileinfo['filepath']
        self.process_file(path)

    def update_grid(self, images):
        self.fv.assert_gui_thread()
        for name, image in images:
            channel = self.fv.get_channel(name)
            viewer = channel.viewer
            canvas = viewer.get_canvas()
            canvas.delete_all_objects()
            #viewer.set_image(image)
            channel.add_image(image)

            self.fv.update_pending()

        if self.tbl_go is not None:
            self.cbar.set_range(self.mag_min, self.mag_max)

            self.plot_stars()

    def set_1k(self, image):
        self.fv.assert_gui_thread()

        channel = self.fv.get_channel_on_demand(self.chname)
        channel.add_image(image)

    def process_file(self, path):
        self.fv.assert_nongui_thread()

        start_time = time.time()
        #opener = PyFitsFileHandler(self.logger)
        opener = self.opener
        try:
            opener.open_file(path)

        except Exception as e:
            self.logger.error("Failed to process image: {}".format(e),
                              exc_info=True)
            return

        if self.current_file is not None:
            self.remove_file(self.current_file)
        self.current_file = path
        self.tbl_go = None
        self.tbl_do = None
        self.tbl_io = None

        images = []
        for dct in opener.hdu_info:
            if dct['name'] == 'PRIMARY':
                # anything we need to do here?  Something with the header?
                pass

            elif dct['name'].startswith('CAM'):
                # load camera image
                idx, name = dct['index'], dct['name']

                image = opener.get_hdu(idx)

                if name == 'CAM1':
                    image.set(path=path)
                    self.fv.gui_do(self.set_1k, image)

                # perform any desired subtractions
                image = self.quick_data_reduce(image, name)
                #image.set(path=f"{path}[{idx}]")

                # update the image in the channel viewer
                images.append((name, image))

            elif dct['name'] == 'detected_objects':
                atbl = opener.get_hdu(dct['index'])
                self.tbl_do = atbl.get_data()

            elif dct['name'] == 'guide_objects':
                atbl = opener.get_hdu(dct['index'])
                self.tbl_go = atbl.get_data()

            elif dct['name'] == 'identified_objects':
                atbl = opener.get_hdu(dct['index'])
                self.tbl_io = atbl.get_data()

            else:
               self.logger.info("Unrecognized HDU: name='{}'".format(dct['name']))

        # determine max and min magnitude
        if self.tbl_go is not None:
            mags = self.tbl_go['mag']
            self.mag_max = np.max(mags)
            self.mag_min = np.min(mags)

        opener.close()
        end_time = time.time()
        self.logger.info("processing time %.4f sec" % (end_time - start_time))

        self.fv.gui_do_oneshot('pfsag_update', self.update_grid, images)

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
        opener = PyFitsFileHandler(self.logger)
        try:
            opener.open_file(path)

        except Exception as e:
            self.logger.error("Failed to open calib file '{}': {}".format(path, e),
                              exc_info=True)
            return

        for dct in opener.hdu_info:
            if dct['name'].startswith('CAM'):
                # load camera image
                idx, name = dct['index'], dct['name']

                image = opener.get_hdu(idx)
                dct[name] = image.get_data()

        opener.close()

    def plot_stars(self):
        if self.tbl_io is None:
            return
        mags = self.tbl_go['mag']

        # delete previously plotted objects
        for cam_num in (1, 2, 3, 4, 5, 6):
            cam_id = 'CAM{}'.format(cam_num)
            channel = self.fv.get_channel(cam_id)
            viewer = channel.viewer
            canvas = viewer.get_canvas()
            canvas.delete_objects_by_tag(canvas.get_tags_by_tag_pfx('_io'))

        # plot identified objects
        for io_idx, io_row in enumerate(self.tbl_io):
            do_row = self.tbl_do[io_row['detected_object_id']]
            cam_num = do_row['camera_id']
            ctr_x, ctr_y = do_row['centroid_x'], do_row['centroid_y']

            # camera indexes are now 0-based, while HDUs are numbered from 1
            cam_id = 'CAM{}'.format(cam_num + 1)
            channel = self.fv.get_channel(cam_id)
            viewer = channel.viewer
            canvas = viewer.get_canvas()

            _go_row_num = int(io_row['guide_object_id'])
            mag = mags[_go_row_num]

            # add circle for detected position
            color = self.get_color(mag)
            radius = 15
            objs = []

            if self.settings.get('plot_guide_stars', False):
                c = self.dc.Circle(ctr_x, ctr_y, radius,
                                   color=color, linewidth=2)
                p = self.dc.Point(ctr_x, ctr_y, radius, style='plus',
                                  color=color, linewidth=2)
                objs.extend([c, p])

            if self.settings.get('plot_offsets', False):
                gde_x, gde_y = (io_row['guide_object_xdet'],
                                io_row['guide_object_ydet'])
                # ra_deg = self.tbl_go['ra'][_go_row_num]
                # dec_deg = self.tbl_go['dec'][_go_row_num]
                # image = viewer.get_image()
                # gde_x, gde_y = image.radectopix(ra_deg, dec_deg)

                c = self.dc.Circle(gde_x, gde_y, radius,
                                   color=color, linestyle='dash',
                                   linewidth=2)
                l = self.dc.Line(ctr_x, ctr_y, gde_x, gde_y,
                                 color=color, linestyle='solid',
                                 linewidth=2, arrow='end')
                objs.extend([c, l])

            canvas.add(self.dc.CompoundObject(*objs),
                       tag=f'_io{io_idx}', redraw=False)

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

    def toggle_plot_guide_stars_cb(self, w, tf):
        self.settings.set(plot_guide_stars=tf)
        self.plot_stars()

    def toggle_plot_offsets_cb(self, w, tf):
        self.settings.set(plot_offsets=tf)
        self.plot_stars()

    def subtract_bias_cb(self, w, tf):
        self.settings.set(subtract_bias=tf)
        if self.current_file is not None:
            self.process_file(self.current_file)

    def subtract_bg_cb(self, w, tf):
        self.settings.set(subtract_background=tf)
        if self.current_file is not None:
            self.process_file(self.current_file)

    def subtract_dark_cb(self, w, tf):
        self.settings.set(subtract_dark=tf)
        if self.current_file is not None:
            self.process_file(self.current_file)

    def divide_flat_cb(self, w, tf):
        self.settings.set(divide_flat=tf)
        if self.current_file is not None:
            self.process_file(self.current_file)

    def set_flats_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.flat)

    def set_darks_cb(self, w):
        path = w.get_text().strip()
        self.read_calib(path, self.dark)

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

                        start_time = time.time()
                        if start_time < self.last_image_time + self.rate_limit:
                            self.logger.info(f"skipping file '{filename}' for rate limit")
                            self.remove_file(filepath)
                            continue
                        self.last_image_time = start_time

                        self.logger.info(f"new file detected: '{filename}'")
                        if self.pause_flag:
                            self.logger.info("plugin is paused, skipping new file")
                            continue
                        fits_tot.append(filename)

                        #self.fv.nongui_do(self.fv.load_file, filepath,
                        #                  chname=self.chname)
                        if filepath in self.processed:
                            self.logger.info(f"file '{filepath}' has already been processed")
                            continue
                        self.processed.add(filepath)
                        self.fv.nongui_do(self.process_file, filepath)
            self.logger.info("---------")
            self.logger.info("{} files: {}".format(len(loop_tot), loop_tot))
            self.logger.info("{} fits: {}".format(len(fits_tot), fits_tot))

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

    def pause_cb(self, w, tf):
        self.pause_flag = tf

    def set_rate_limit_cb(self, w, val):
        self.rate_limit = val

    def __str__(self):
        return 'PFS_AG'
