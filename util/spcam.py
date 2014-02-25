#
# spcam.py -- Suprime-Cam data processing routines 
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import numpy

#from ginga import AstroImage
from ginga.misc import Bunch


def get_regions(image):
    """Extract the keywords defining the overscan and effective pixel
    regions in a SPCAM image.  The data is returned in a dictionary of
    bunches.  The keys of the dictionary are the channel numbers, plus
    'image'.
    """
    wd, ht = image.get_size()
    d = {}
    xcut = 0
    newwd = 0
    l = []
    for channel in (1, 2, 3, 4):
        base = 'S_EF'
        efminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
        efmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
        efminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
        efmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
        base = 'S_OS'
        osminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
        osmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
        osminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
        osmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
        xcut += osmaxx - osminx + 1
        newwd += efmaxx + 1 - efminx 
        d[channel] = Bunch.Bunch(
            efminx=efminx, efmaxx=efmaxx, efminy=efminy, efmaxy=efmaxy,
            osminx=osminx, osmaxx=osmaxx, osminy=osminy, osmaxy=osmaxy)
        l.append(d[channel])

    # figure out starting x position of channel within image
    l.sort(cmp=lambda x, y: x.efmaxx - y.efmaxx)
    startposx = 0
    for ch in l:
        ch.setvals(startposx=startposx)
        startposx += ch.efmaxx + 1 - ch.efminx
    
    ycut = osmaxy - osminy + 1
    newht = efmaxy + 1 - efminy 
    d['image']=Bunch.Bunch(xcut=xcut, ycut=ycut,
                           newwd=newwd, newht=newht)

    return d


def subtract_overscan(data_np, d):
    """Remove the overscan regions from a SPCAM image data array (data_np).
    'd' is a dictionary of information about the overscan and effective
    pixel regions as returned by get_regions().

    Output:
      A new, smaller array is returned.
    """

    # create new output array the size of the sum of the image
    # effective pixels
    info = d['image']
    newwd, newht = info.newwd, info.newht
    print "effective pixel size %dx%d" % (newwd, newht)
    out = numpy.empty((newht, newwd), dtype=float)

    # original image size
    ht, wd = data_np.shape[:2]

    for channel in (1, 2, 3, 4):
        print "processing channel %d" % (channel)
        ch = d[channel]

        # get median of each row in overscan area for this channel
        ovsc_median = numpy.median(data_np[ch.efminy:ch.efmaxy+1,
                                       ch.osminx:ch.osmaxx+1], axis=1)
        # calculate size of effective pixels area for this channel
        efwd = ch.efmaxx + 1 - ch.efminx
        efht = ch.efmaxy + 1 - ch.efminy
        len_ovsc = ovsc_median.shape[0]

        assert len_ovsc == efht, \
               ValueError("median array len (%d) doesn't match effective pixel len (%d)" % (
            len_ovsc, efht))

        ovsc_median = ovsc_median.reshape((efht, 1))
        ovsc_median = numpy.repeat(ovsc_median, efwd, axis=1)

        j = ch.startposx

        # Cut effective pixel region into output array
        xlo, xhi, ylo, yhi = j, j + efwd, 0, efht
        out[ylo:yhi, xlo:xhi] = data_np[ch.efminy:ch.efmaxy+1,
                                        ch.efminx:ch.efmaxx+1]
        # Subtract overscan medians
        #print ovsc_median
        #out[ylo:yhi, xlo:xhi] -= ovsc_median[:, None]
        out[ylo:yhi, xlo:xhi] -= ovsc_median

    return out


def step_two(image):

    d = get_regions(image)
    try:
        newarr = subtract_overscan(image.get_data(), d)
    except Exception as e:
        print "Error: %s" % (str(e))
    image.set_data(newarr)
    return image


#END
