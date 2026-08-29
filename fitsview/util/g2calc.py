#
# g2calc.py -- Gen2 image quality calculations on FITS data
#
# E. Jeschke
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#

import time


from ginga.misc import Bunch
from ginga.util import iqcalc

# special for fitsview & guideview.  "eclipse" is a compiled extension that
# is often not built; its absence must not stop this module from importing,
# since every plugin that measures a star imports it.
try:
    from eclipse import qualsize as eclipse_qualsize
except ImportError:
    eclipse_qualsize = None

# Star selection limits compiled into the SOSS qualsize(), from qualsize.c.
# The fallback uses them so that it accepts the same candidates.
MINFWHM = 0.1
MAXFWHM = 50.0
MINELIPS = 0.0
EDGEW = 0.01


class IQCalc(iqcalc.IQCalc):

    def qualsize_old(self, image, x1=None, y1=None, x2=None, y2=None,
                     radius=5, threshold=None,
                     minfwhm=MINFWHM, maxfwhm=MAXFWHM,
                     minelipse=MINELIPS, edgew=EDGEW):
        """Measure a star with Kosugi-san's original SOSS algorithm.

        Falls back to :meth:`ginga.util.iqcalc.IQCalc.qualsize` where the
        "eclipse" extension providing that algorithm is not built.  The
        result carries the same ``x``, ``y``, ``fwhm``, ``brightness``,
        ``objx`` and ``objy`` either way, so callers need not care which
        one ran; the values themselves do differ, as the two algorithms
        measure stars differently.

        ``minfwhm``, ``maxfwhm``, ``minelipse`` and ``edgew`` select which
        candidates are acceptable.  They default to the values the SOSS
        algorithm has compiled in, so that the fallback accepts the same
        stars it would have, rather than IQCalc's stricter defaults --
        those reject small and elongated objects the SOSS routine measures
        happily.  The SOSS routine cannot be told anything else, so these
        only take effect on the fallback.
        """
        if eclipse_qualsize is None:
            self.logger.debug("eclipse is not available; measuring with "
                              "IQCalc.qualsize() instead")
            return self.qualsize(image, x1=x1, y1=y1, x2=x2, y2=y2,
                                 radius=radius, threshold=threshold,
                                 minfwhm=minfwhm, maxfwhm=maxfwhm,
                                 minelipse=minelipse, edgew=edgew)

        x1, y1, x2, y2 = int(x1), int(y1), int(x2), int(y2)
        data = image.cutout_data(x1, y1, x2, y2, astype='float32')

        start_time = time.time()
        (x, y, fwhm, brightness, skylevel,
         objx, objy) = eclipse_qualsize.qualsize(data)
        # NOTE: fwhm_x and fwhm_y carry the one combined value the SOSS
        # routine measures -- it does not resolve the two axes separately.
        # They are here so that the result has the same shape whichever
        # algorithm ran, since callers pair them with CDELT1 and CDELT2.
        # For non-square pixels this can only give the mean of the two
        # scales; use the IQCalc algorithm to measure the axes properly.
        qs = Bunch.Bunch(x=x, y=y, fwhm=fwhm, fwhm_x=fwhm, fwhm_y=fwhm,
                         brightness=brightness,
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
