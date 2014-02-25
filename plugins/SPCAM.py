#
# SPCAM.py -- Suprime-Cam plugin for Ginga FITS viewer
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import Queue

from astro.frame import Frame

from Gen2.fitsview.util import spcam
from ginga.misc.plugins import Mosaic


class SPCAM(Mosaic.Mosaic):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(SPCAM, self).__init__(fv, fitsimage)

        # Set preferences for destination channel
        prefs = self.fv.get_preferences()
        self.settings = prefs.createCategory('plugin_SPCAM')
        self.settings.setDefaults(annotate_images=False, fov_deg=0.72,
                                  match_bg=False, trim_px=0,
                                  merge=True, num_threads=4)
        self.settings.load(onError='silent')

        self.queue = Queue.Queue()
        self.current_exp_num = 0
        
        # define mosaic preprocessing step for SPCAM 
        self.set_preprocess(spcam.step_two)

        self.fv.set_callback('file-notify', self.file_notify_cb)

        self.timer = self.fv.get_timer()
        self.timer.add_callback('expired', self.process_frames)
        self.process_interval = 0.5

        self.mosaic_chname = 'SPCAM_Online'


    def instructions(self):
        self.tw.set_text("""Frames will be mosaiced as they arrive.""")
            
    def get_latest_frames(self, pathlist):
        new_frlist = []
        new_exposure = False

        for path in pathlist:
            frame = Frame(path=path)
            # if not a SPCAM frame then drop it
            if frame.inscode != 'SUP':
                continue

            # if frame number doesn't belong to current exposure
            # then drop it
            if frame.number < self.current_exp_num:
                continue

            exp_num = (frame.number // 10) * 10
            if exp_num > self.current_exp_num:
                # There is a new exposure
                self.current_exp_num = exp_num
                new_frlist = [ path ]
                new_exposure = True
            else:
                new_frlist.append(path)

        return (new_frlist, new_exposure)
    
            
    def process_frames(self, timer):
        self.logger.debug("processing queued frames")

        # Get all files stored in the queue
        paths = []
        while True:
            try:
                path = self.queue.get(block=False)
                paths.append(path)
            except Queue.Empty:
                break

        if len(paths) == 0:
            return

        if self.gui_up:
            paths, new_mosaic = self.get_latest_frames(paths)

            self.mosaic(paths, new_mosaic=new_mosaic)

            
    def file_notify_cb(self, fv, path):
        self.logger.debug("file notify: %s" % (path))
        self.queue.put(path)
        self.timer.set(self.process_interval)
        
        
    def __str__(self):
        return 'spcam'
    

#END
