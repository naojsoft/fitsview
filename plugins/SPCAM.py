#
# SPCAM.py -- Suprime-Cam plugin for Ginga FITS viewer
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import os, re, glob
import threading, Queue

from ginga import AstroImage
from ginga.misc.plugins import Mosaic
from ginga.misc import Widgets
from ginga.util import dp

from astro.frame import Frame
from Gen2.fitsview.util import spcam


class SPCAM(Mosaic.Mosaic):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(SPCAM, self).__init__(fv, fitsimage)

        # Set preferences for destination channel
        prefs = self.fv.get_preferences()
        self.settings = prefs.createCategory('plugin_SPCAM')
        self.settings.setDefaults(annotate_images=False, fov_deg=0.72,
                                  match_bg=False, trim_px=0,
                                  merge=True, num_threads=4,
                                  drop_creates_new_mosaic=True,
                                  use_flats=False, flat_dir='')
        self.settings.load(onError='silent')

        self.queue = Queue.Queue()
        self.current_exp_num = 0
        
        # define mosaic preprocessing step for SPCAM 
        self.set_preprocess(self.mangle_image)

        self.fv.enable_callback('file-notify')
        self.fv.add_callback('file-notify', self.file_notify_cb)

        self.timer = self.fv.get_timer()
        self.timer.add_callback('expired', self.process_frames)
        self.process_interval = 0.2
        
        self.dr = spcam.SuprimeCamDR(logger=self.logger)
        
        #self.mosaic_chname = 'SPCAM_Online'
        # For flat fielding
        self.flat = {}


    def build_gui(self, container):
        super(SPCAM, self).build_gui(container)

        vbox = self.w.vbox
        
        fr = Widgets.Frame("Flats")

        captions = [
            ("Use flats", 'checkbutton'),
            ("Flat dir:", 'label', 'flat_dir', 'entry'),
            ("Load Flats", 'button'),
            ]
        w, b = Widgets.build_info(captions)
        self.w.update(b)

        b.flat_dir.set_length(32)
        b.flat_dir.set_text(self.settings.get('flat_dir', ''))
        b.load_flats.add_callback('activated', self.load_flats_cb)
        b.use_flats.set_tooltip("Flat field tiles as they arrive")
        use_flats = self.settings.get('use_flats', False)
        b.use_flats.set_state(use_flats)
        b.flat_dir.set_tooltip("Directory containing flat field tiles")
        b.load_flats.set_tooltip("Load flat field tiles from directory")

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)


    def instructions(self):
        self.tw.set_text("""Frames will be mosaiced as they arrive.""")
            
    def get_latest_frames(self, pathlist):
        new_frlist = []
        new_exposure = False

        for path in pathlist:
            path = self.fv.get_filepath(path)
            #print "path is", path
                
            frame = Frame(path=path)
            # if not an instrument frame then drop it
            if frame.inscode != self.dr.inscode:
                continue

            # if frame number doesn't belong to current exposure
            # then drop it
            if frame.number < self.current_exp_num:
                continue

            exp_num = (frame.number // self.dr.num_frames) * self.dr.num_frames
            if exp_num > self.current_exp_num:
                # There is a new exposure
                self.current_exp_num = exp_num
                new_frlist = [ path ]
                new_exposure = True
            else:
                new_frlist.append(path)

        return (new_frlist, new_exposure)
    
            
    def process_frames(self, timer):
        self.logger.info("processing queued frames")

        # Get all files stored in the queue
        paths = []
        while True:
            try:
                path = self.queue.get(block=False)
                paths.append(path)
            except Queue.Empty:
                break

        self.logger.info("debug=%s" % str(paths))
        
        if len(paths) == 0:
            return

        if self.gui_up:
            paths, new_mosaic = self.get_latest_frames(paths)

            self.mosaic(paths, new_mosaic=new_mosaic)

    def file_notify_cb(self, fv, path):
        self.logger.debug("file notify: %s" % (path))
        self.queue.put(path)
        self.timer.cond_set(self.process_interval)

    ## def drop_cb(self, canvas, paths):
    ##     self.logger.info("files dropped: %s" % str(paths))
    ##     if self.gui_up:
    ##         paths, new_mosaic = self.get_latest_frames(paths)

    ##         self.fv.nongui_do(self.fv.error_wrap, self.mosaic, paths,
    ##                           new_mosaic=new_mosaic)
    ##     return True
        
    def mangle_image(self, image):
        d = self.dr.get_regions(image)
        header = {}
        data_np = image.get_data()

        # subtract overscan region
        result = self.dr.subtract_overscan_np(data_np, d,
                                              header=header)

        # flat field this piece, if flat provided
        do_flat = self.w.use_flats.get_state()
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

    def __str__(self):
        return 'spcam'
    

#END
