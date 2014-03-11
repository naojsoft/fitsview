#
# HSC.py -- Suprime-Cam plugin for Ginga FITS viewer
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
from ginga.misc import Widgets
from ginga.util import dp

from Gen2.fitsview.util import hsc
from Gen2.fitsview.plugins import SPCAM


class HSC(SPCAM.SPCAM):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(HSC, self).__init__(fv, fitsimage)

        # Set preferences for destination channel
        prefs = self.fv.get_preferences()
        self.settings = prefs.createCategory('plugin_HSC')
        self.settings.setDefaults(annotate_images=False, fov_deg=2.0,
                                  match_bg=False, trim_px=0,
                                  merge=True, num_threads=6,
                                  drop_creates_new_mosaic=True)
        self.settings.load(onError='silent')

        self.dr = hsc.HyperSuprimeCamDR(logger=self.logger)
        
        #self.mosaic_chname = 'HSC_Online'


    def __str__(self):
        return 'hsc'
    

#END
