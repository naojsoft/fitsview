#
# SPCAM.py -- Suprime-Cam quick look plugin for Ginga FITS viewer
#
# E. Jeschke
#
#
"""
A plugin for the Ginga scientific image viewer for quick look viewing and
mosaicing Suprime-Cam data.

Installation:
  $ mkdir $HOME/.ginga
  $ mkdir $HOME/.ginga/plugins
  $ cp SPCAM.py $HOME/.ginga/plugins

Running:
  $ ginga [other options] --plugins=SPCAM ...

NOTE: this requires the "naojutils" module, available at
          https://github.com/naojsoft/naojutils
"""
import os, re, glob
import time
import threading
import queue as Queue

from ginga import AstroImage
from ginga.rv.plugins import Mosaic
from ginga.misc import Bunch, Future
from ginga.util import dp, iohelper
from ginga.gw import Widgets

# You should install module "naojutils" to run this plugin.
#   Get it here--> https://github.com/naojsoft/naojutils
from naoj.frame import Frame
from naoj.spcam import spcam_dr


class SPCAM(Mosaic.Mosaic):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(SPCAM, self).__init__(fv, fitsimage)

        # Set preferences for destination channel
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_SPCAM')
        self.settings.set_defaults(annotate_images=False, fov_deg=0.72,
                                   match_bg=False, trim_px=0,
                                   merge=True, num_threads=4,
                                   drop_creates_new_mosaic=False,
                                   mosaic_hdus=False, skew_limit=0.1,
                                   allow_expand=False, expand_pad_deg=0.01,
                                   use_flats=False, flat_dir='',
                                   mosaic_new=True, make_thumbs=True,
                                   reuse_image=False)
        self.settings.load(onError='silent')

        self.queue = Queue.Queue()
        self.current_exp_num = 0
        self.sub_bias = False

        # map of exposures -> list of paths
        self.exp_pathmap = {}

        # define mosaic preprocessing step for SPCAM
        self.set_preprocess(self.mangle_image)

        self.fv.enable_callback('file-notify')
        self.fv.add_callback('file-notify', self.file_notify_cb)

        self.timer = self.fv.get_timer()
        self.timer.add_callback('expired', self.process_frames)
        self.process_interval = 2.0

        self.dr = spcam_dr.SuprimeCamDR(logger=self.logger)
        self.basedir = "/gen2/share/data/SPCAM"

        # mosaic results channel
        self.mosaic_chname = 'SPCAM_Online'

        # For flat fielding
        self.flat = {}
        self.use_flats = self.settings.get('use_flats', False)


    def build_gui(self, container):
        super(SPCAM, self).build_gui(container)

        if not self.fv.has_channel(self.mosaic_chname):
            self.fv.add_channel(self.mosaic_chname)

        vbox = self.w.vbox

        fr = Widgets.Frame("Bias")

        captions = [
            ("Sub Bias", 'checkbutton'),
            ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.sub_bias.set_tooltip("Subtract bias calculated from overscan regions")
        b.sub_bias.set_state(self.sub_bias)
        b.sub_bias.add_callback('activated', self.sub_bias_cb)

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        fr = Widgets.Frame("Flats")

        captions = [
            ("Use flats", 'checkbutton'),
            ("Flat dir:", 'label', 'flat_dir', 'entry'),
            ("Load Flats", 'button'),
            ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        #b.flat_dir.set_length(512)
        b.flat_dir.set_text(self.settings.get('flat_dir', ''))
        b.load_flats.add_callback('activated', self.load_flats_cb)
        b.use_flats.set_tooltip("Flat field tiles as they arrive")
        b.use_flats.set_state(self.use_flats)
        b.use_flats.add_callback('activated', self.use_flats_cb)
        b.flat_dir.set_tooltip("Directory containing flat field tiles")
        b.load_flats.set_tooltip("Load flat field tiles from directory")

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        fr = Widgets.Frame("Load Exposure")

        captions = [
            ("Load exp:", 'label', 'exposure', 'entry',
             "Load Exp", 'button'),
            ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.load_exp.add_callback('activated', self.load_exp_cb)
        b.exposure.add_callback('activated', self.load_exp_cb)

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        fr = Widgets.Frame("Quick Color Maps")

        captions = [
            ("Gray", 'button', "JT", 'button'),
            ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.gray.add_callback('activated', lambda w: self.set_cmap('gray'))
        b.jt.add_callback('activated', lambda w: self.set_cmap('jt'))

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)


    def instructions(self):
        self.tw.set_text("""Frames will be mosaiced as they arrive.""")

    def get_exp_num(self, frame):
        exp_num = (frame.number // self.dr.num_frames) * self.dr.num_frames
        return exp_num

    def get_latest_frames(self, pathlist):
        new_frlist = []
        new_exposure = False
        imname = None
        exposures = set([])

        for path in pathlist:
            info = iohelper.get_fileinfo(path)
            self.logger.info("getting path")
            path = info.filepath
            self.logger.info("path is %s" % (path))

            frame = Frame(path=path)
            # if not an instrument frame then drop it
            if frame.inscode != self.dr.inscode:
                continue

            # calculate exposure number
            #exp_num = (frame.number // self.dr.num_frames) * self.dr.num_frames
            exp_num = self.get_exp_num(frame)

            # add paths to exposure->paths map
            exp_id = Frame(path=path)
            exp_id.number = exp_num
            exp_frid = str(exp_id)
            bnch = Bunch.Bunch(paths=set([]), name=exp_frid)
            exp_bnch = self.exp_pathmap.setdefault(exp_frid, bnch)
            exp_bnch.paths.add(path)
            if len(exp_bnch.paths) == 1:
                exp_bnch.setvals(typical=path, added_to_contents=False)
            exposures.add(exp_frid)

            # if frame number doesn't belong to current exposure
            # then drop it
            if frame.number < self.current_exp_num:
                continue

            if exp_num > self.current_exp_num:
                # There is a new exposure
                self.current_exp_num = exp_num
                new_frlist = [ path ]
                new_exposure = True
                imname = exp_frid
            else:
                new_frlist.append(path)

        return (new_frlist, new_exposure, imname, exposures)


    def mk_future(self, bnch):
        def load_mosaic(bnch):
            self.fv.assert_nongui_thread()
            paths = list(bnch.paths)
            return self.mosaic(paths,
                               new_mosaic=True, name=bnch.name)
        future = Future.Future()
        future.freeze(load_mosaic, bnch)
        return future

    def add_to_channel(self, exposures):
        self.logger.info("adding to channel: exposures=%s" % (str(exposures)))
        try:
            channel = self.fv.get_channel(self.mosaic_chname)

            for exp_frid in exposures:
                # get the information about this exposure
                exp_bnch = self.exp_pathmap[exp_frid]
                if exp_bnch.added_to_contents:
                    continue

                self.logger.debug("Exposure '%s' not yet added to channel" % (
                    exp_frid))

                # load the representative image
                #image = AstroImage.AstroImage(logger=self.logger)
                # TODO: is this load even necessary?  Would be good
                # if we could just load the headers
                path = exp_bnch.typical

                # make a future that will load the mosaic if we need
                # to recreate it
                image_future = self.mk_future(exp_bnch)

                info = Bunch.Bunch(name=exp_frid, path=None,
                                   image_future=image_future)

                exp_bnch.added_to_contents = True

                # add this to the channel contents
                self.logger.info("adding to channel: info=%s" % (str(info)))
                channel.add_image_info(info)

        except Exception as e:
            self.logger.warn("Error adding to channel: %s" % (str(e)))
            return

    def process_frames(self, timer):
        self.fv.gui_do(self._process_frames)

    def _process_frames(self):
        self.fv.assert_gui_thread()
        self.logger.info("processing queued frames")

        # Get all files stored in the queue
        paths = []
        while True:
            try:
                path = self.queue.get(block=False)
                paths.append(path)
            except Queue.Empty:
                break

        self.logger.debug("1. paths=%s" % str(paths))
        if len(paths) == 0:
            return

        self.logger.info("paths are: %s" % (str(paths)))
        try:
            paths, new_mosaic, imname, exposures = self.get_latest_frames(paths)

            if len(exposures) > 0:
                self.add_to_channel(exposures)

            self.logger.debug("2. paths=%s" % str(paths))
            if len(paths) == 0:
                return

        except Exception as e:
            self.logger.error("error adding to channel: %s" % (str(e)))

        mosaic_new = self.settings.get('mosaic_new', False)
        self.logger.info("mosaic_new=%s new_mosaic=%s imname=%s" % (
            mosaic_new, new_mosaic, imname))
        if self.gui_up and mosaic_new:
            self.logger.info("mosaicing %s" % (str(paths)))
            self.fv.nongui_do(self.mosaic, paths, name=imname,
                              new_mosaic=new_mosaic)

    def file_notify_cb(self, fv, path):
        self.logger.debug("file notify: %s" % (path))
        self.queue.put(path)
        self.timer.cond_set(self.process_interval)

    def drop_cb(self, canvas, paths):
        self.logger.info("files dropped: %s" % str(paths))
        for path in paths:
            self.queue.put(path)
        self.timer.cond_set(self.process_interval)
        return True

    def load_exp_cb(self, w):
        frid = self.w.exposure.get_text().strip()
        path = os.path.join(self.basedir, frid)
        paths = self.dr.get_file_list(path)
        self.logger.info("paths are: %s" % (paths))

        new_paths, new_mosaic, imname, exposures = self.get_latest_frames(paths)
        if len(exposures) > 0:
            self.add_to_channel(exposures)

        self.fv.nongui_do(self.mosaic, paths, name=imname, new_mosaic=True)

    def set_cmap(self, cmname):
        self.fitsimage.set_color_map(cmname)

    def mangle_image(self, image):
        d = self.dr.get_regions(image)
        header = {}
        data_np = image.get_data()

        # subtract overscan region
        result = self.dr.subtract_overscan_np(data_np, d,
                                              header=header,
                                              sub_bias=self.sub_bias)

        # flat field this piece, if flat provided
        do_flat = self.use_flats
        if do_flat and (len(self.flat) > 0):
            try:
                ccd_id = int(image.get_keyword('DET-ID'))
                result /= self.flat[ccd_id]
            except Exception as e:
                self.logger.warn("Error applying flat field: %s" % (str(e)))

        newimage = dp.make_image(result, image, header)
        return newimage

    def _load_flats(self, datadir):
        # TODO: parallelize this
        self.fv.assert_nongui_thread()

        path_glob = os.path.join(datadir, '*-*.fits')
        d = {}
        paths = glob.glob(path_glob)
        if len(paths) != self.dr.num_ccds:
            self.fv.gui_do(self.fv.show_error, "Number of flat files (%d) does not match number of CCDs (%d)" % (
                len(paths), self.dr.num_frames))
            return

        self.update_status("Loading flats...")
        self.init_progress()

        self.total_files = max(1, len(paths))
        self.ingest_count = 0

        num_threads = self.settings.get('num_threads', 4)
        groups = self.split_n(paths, num_threads)
        for group in groups:
            self.fv.nongui_do(self._load_some, d, group)

    def _load_some(self, d, paths):
        for path in paths:
            match = re.match(r'^.+\-(\d+)\.fits$', path)
            if match:
                ccd_id = int(match.group(1))
                image = AstroImage.AstroImage(logger=self.logger)
                image.load_file(path)

                with self.lock:
                    d[ccd_id] = image.get_data()
                    self.ingest_count += 1
                    count = self.ingest_count

                self.update_progress(float(count)/self.total_files)

        if count == self.total_files:
            self.flat = d
            self.end_progress()
            self.update_status("Flats loaded.")

    def load_flats_cb(self, w):
        dirpath = self.w.flat_dir.get_text().strip()
        # Save the setting
        self.settings.set(flat_dir=dirpath)

        self.fv.nongui_do(self._load_flats, dirpath)

    def use_flats_cb(self, w, tf):
        self.use_flats = tf

    def sub_bias_cb(self, w, tf):
        self.sub_bias = tf

    def __str__(self):
        return 'spcam'


#END
