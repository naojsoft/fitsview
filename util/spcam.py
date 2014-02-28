#
# spcam.py -- Suprime-Cam data processing routines 
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import os.path
import numpy

from ginga.misc import Bunch
from ginga.util import dp


def get_regions(image, pfx='S'):
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
        base = pfx + '_EF'
        efminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
        efmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
        efminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
        efmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
        base = pfx + '_OS'
        osminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
        osmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
        osminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
        osmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
        xcut += osmaxx - osminx + 1
        newwd += efmaxx + 1 - efminx 

        gain = float(image.get_keyword("%s_GAIN%d" % (pfx, channel)))
        d[channel] = Bunch.Bunch(
            efminx=efminx, efmaxx=efmaxx, efminy=efminy, efmaxy=efmaxy,
            osminx=osminx, osmaxx=osmaxx, osminy=osminy, osmaxy=osmaxy,
            gain=gain)
        l.append(d[channel])

    # figure out starting x position of channel within image
    l.sort(cmp=lambda x, y: x.efmaxx - y.efmaxx)
    startposx = 0
    for ch in l:
        ch.setvals(startposx=startposx)
        startposx += ch.efmaxx + 1 - ch.efminx
    
    ycut = osmaxy - osminy + 1
    newht = efmaxy + 1 - efminy 
    d['image'] = Bunch.Bunch(xcut=xcut, ycut=ycut,
                             newwd=newwd, newht=newht)

    return d


def subtract_overscan_np(data_np, d, pfx='S', header=None):
    """Subtract the median bias calculated from the overscan regions
    from a SPCAM image data array.  The resulting image is trimmed to
    remove the overscan regions.

    Parameters
    ----------
    data_np: numpy array
        a 2D data array of pixel values
    d: dict
        a dictionary of information about the overscan and effective
        pixel regions as returned by get_regions().

    Returns:
    out: numpy array
        a new, smaller array with the result data
    """

    # create new output array the size of the sum of the image
    # effective pixels
    info = d['image']
    newwd, newht = info.newwd, info.newht
    #print "effective pixel size %dx%d" % (newwd, newht)
    out = numpy.empty((newht, newwd), dtype=float)
    if header != None:
        header['NAXIS1'] = newwd
        header['NAXIS2'] = newht

    # original image size
    ht, wd = data_np.shape[:2]

    for channel in (1, 2, 3, 4):
        #print "processing channel %d" % (channel)
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
        #out[ylo:yhi, xlo:xhi] -= ovsc_median[:, None]
        out[ylo:yhi, xlo:xhi] -= ovsc_median

        # Update header for effective regions
        if header != None:
            base = pfx + '_EF'
            header["%sMN%d1" % (base, channel)] = xlo + 1
            header["%sMX%d1" % (base, channel)] = xhi + 1
            header["%sMN%d2" % (base, channel)] = ylo + 1
            header["%sMX%d2" % (base, channel)] = yhi + 1

    return out


## def make_flat_v1(flatlist, bias=None, pfx='S'):

##     flats = []
##     for path in flatlist:
##         image = AstroImage.AstroImage()
##         image.load_file(path)

##         data_np = image.get_data()
##         # TODO: subtract optional bias image

##         # subtract overscan and trim
##         d = get_regions(image, pfx=pfx)
##         header = {}
##         newarr = subtract_overscan_np(data_np, d, header=header)
        
##         flats.append(newarr)

##     # Take the median of the individual frames
##     flat = numpy.median(numpy.array(flats), axis=0)
##     #print flat.shape

##     # Normalize flat
##     flat = flat / numpy.mean(flat.flat)

##     img_flat = dp.make_image(flat, image, header)
##     return img_flat


## def make_flats(datadir, explist, pfx='S'):

##     flats = []
##     for i in xrange(10):
##         flatlist = []
##         for exp in explist:
##             path = os.path.join(datadir, exp.upper()+'.fits')
##             frame = Frame(path=path)
##             frame.number += i
##             flatlist.append(os.path.join(datadir, str(frame)+'.fits'))

##         flats.append(make_flat(flatlist, pfx=pfx))

##     return flats
            
    
def step2(image, pfx='S'):

    d = get_regions(image, pfx=pfx)
    header = {}
    data_np = image.get_data()

    result = subtract_overscan_np(data_np, d, header=header)

    newimage = dp.make_image(result, image, header)
    return newimage


def step3(image, flat):

    data_np = image.get_data()
    flat_np = flat.get_data()

    result = data_np / flat_np

    newimage = dp.make_image(result, image, {})
    return newimage


#END
