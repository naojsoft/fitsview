#
# g2calc.py -- Gen2 image quality calculations on FITS data
#
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c) Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#

import math
import time

import numpy
    
from ginga.misc import Bunch
from ginga.util import iqcalc

# special for fitsview & guideview
import qualsize


class IQCalc(iqcalc.IQCalc):

    def qualsize_old(self, image, x1=None, y1=None, x2=None, y2=None,
                     radius=5, threshold=None):

        x1, y1, x2, y2 = int(x1), int(y1), int(x2), int(y2)
        data = image.cutout_data(x1, y1, x2, y2, astype='float32')

        start_time = time.time()
        (x, y, fwhm, brightness, skylevel, objx, objy) = qualsize.qualsize(data)
        qs = Bunch.Bunch(x=x, y=y, fwhm=fwhm, brightness=brightness,
                         skylevel=skylevel, objx=objx, objy=objy)
        elapsed = time.time() - start_time
        
        # Add back in offsets into image to get correct values with respect
        # to the entire image
        qs.x += x1
        qs.y += y1
        qs.objx += x1
        qs.objy += y1
        self.logger.debug("obj=%f,%f fwhm=%f sky=%f bright=%f (%f sec)" % (
            qs.objx, qs.objy, qs.fwhm, qs.skylevel, qs.brightness, elapsed))

        return qs
     
#END
