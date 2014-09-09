#
# spcam.py -- Suprime-Cam data processing routines 
# 
# Eric Jeschke (eric@naoj.org)
#
# Copyright (c)  Eric R. Jeschke.  All rights reserved.
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import os
import re, glob

import numpy

from ginga import AstroImage
from ginga.misc import Bunch, log
from ginga.util import dp, mosaic, io_fits

from astro.frame import Frame


class SuprimeCamDR(object):

    def __init__(self, logger=None):
        super(SuprimeCamDR, self).__init__()
        
        if logger == None:
            logger = log.get_logger(level=20, log_stderr=True)
        self.logger = logger

        self.pfx = 'S'
        self.num_ccds = 10
        self.num_frames = 10
        self.inscode = 'SUP'
        self.fov = 0.72
        
    def get_regions(self, image):
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
            base = self.pfx + '_EF'
            efminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
            efmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
            efminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
            efmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
            base = self.pfx + '_OS'
            osminx = int(image.get_keyword("%sMN%d1" % (base, channel))) - 1
            osmaxx = int(image.get_keyword("%sMX%d1" % (base, channel))) - 1
            osminy = int(image.get_keyword("%sMN%d2" % (base, channel))) - 1
            osmaxy = int(image.get_keyword("%sMX%d2" % (base, channel))) - 1
            xcut += osmaxx - osminx + 1
            newwd += efmaxx + 1 - efminx 

            gain = float(image.get_keyword("%s_GAIN%d" % (self.pfx, channel)))
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


    def subtract_overscan_np(self, data_np, d, header=None):
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
                base = self.pfx + '_EF'
                header["%sMN%d1" % (base, channel)] = xlo + 1
                header["%sMX%d1" % (base, channel)] = xhi + 1
                header["%sMN%d2" % (base, channel)] = ylo + 1
                header["%sMX%d2" % (base, channel)] = yhi + 1

        return out


    def make_flat(self, flatlist, bias=None, flat_norm=None,
                  logger=None):

        flats = []
        for path in flatlist:
            image = AstroImage.AstroImage(logger=logger)
            image.load_file(path)

            data_np = image.get_data()
            # TODO: subtract optional bias image

            # subtract overscan and trim
            d = self.get_regions(image)
            header = {}
            newarr = self.subtract_overscan_np(data_np, d,
                                               header=header)

            flats.append(newarr)

        # Take the median of the individual frames
        flat = numpy.median(numpy.array(flats), axis=0)
        #print flat.shape

        # Normalize flat, if normalization term provided
        if flat_norm != None:
            flat = flat / flat_norm

        img_flat = dp.make_image(flat, image, header)
        return img_flat


    def make_flat_tiles(self, datadir, explist, output_pfx='flat',
                        output_dir=None):

        # Get the median values for each CCD image
        flats = []
        for i in xrange(self.num_frames):
            flatlist = []
            for exp in explist:
                path = os.path.join(datadir, exp.upper()+'.fits')
                if not os.path.exists(path):
                    continue
                frame = Frame(path=path)
                frame.number += i
                path = os.path.join(datadir, str(frame)+'.fits')
                if not os.path.exists(path):
                    continue
                flatlist.append(path)

            if len(flatlist) > 0:
                flats.append(self.make_flat(flatlist))

        # Normalize the flats
        # TODO: can we use a running median to speed this up without
        # losing much precision?
        # flatarr = numpy.array([ image.get_data() for image in flats ])
        # mval = numpy.median(flatarr.flat)
        flatarr = numpy.array([ numpy.median(image.get_data())
                                for image in flats ])
        mval = numpy.mean(flatarr)

        d = {}
        for image in flats:
            flat = image.get_data()
            flat /= mval
            # no zero divisors
            flat[flat == 0.0] = 1.0
            ccd_id = int(image.get_keyword('DET-ID'))

            if output_dir == None:
                d[ccd_id] = image
            else:
                # write the output file
                name = '%s-%d.fits' % (output_pfx, ccd_id)
                outfile = os.path.join(output_dir, name)
                d[ccd_id] = outfile
                self.logger.debug("Writing output file: %s" % (outfile))
                try:
                    os.remove(outfile)
                except OSError:
                    pass
                image.save_as_file(outfile)
                
        return d


    def get_flat_name(self, pfx, image):
        hdr = image.get_header()
        kwds = dict([ (kwd, hdr[kwd]) for kwd in ('OBJECT', 'FILTER01',
                                                  'DATE-OBS', 'UT-STR') ])
        match = re.match(r'^(\d\d):(\d\d):(\d\d)\.\d+$', kwds['UT-STR'])
        ut = ''.join(match.groups())
        match = re.match(r'^(\d\d\d\d)\-(\d+)\-(\d+)$', kwds['DATE-OBS'])
        date = ''.join(match.groups())
        fname = '%s-flat-%s-%s-%s-%s.fits' % (pfx, date, ut,
                                              kwds['OBJECT'],
                                              kwds['FILTER01'])
        return fname, kwds


    def make_flat_tiles_exp(self, datadir, expstart, num_exp,
                            output_pfx='flat', output_dir=None):

        path = os.path.join(datadir, expstart.upper()+'.fits')

        # make a list of all the exposure ids
        explist = []
        for i in xrange(num_exp):
            frame = Frame(path=path)
            frame.number += i * self.num_frames
            explist.append(str(frame))

        d = self.make_flat_tiles(datadir, explist,
                                 output_pfx=output_pfx,
                                 output_dir=output_dir)
        return d


    def load_flat_tiles(self, datadir):

        path_glob = os.path.join(datadir, '*-*.fits')
        d = {}
        for path in glob.glob(path_glob):
            match = re.match(r'^.+\-(\d+)\.fits$', path)
            if match:
                ccd_id = int(match.group(1))
                image = AstroImage.AstroImage(logger=self.logger)
                image.load_file(path)

                d[ccd_id] = image.get_data()

        return d

def step2(image):

    dr = SuprimeCamDR()
    d = dr.get_regions(image)
    header = {}
    data_np = image.get_data()

    result = dr.subtract_overscan_np(data_np, d, header=header)

    newimage = dp.make_image(result, image, header)
    return newimage


def step3(image, flat):

    data_np = image.get_data()
    flat_np = flat.get_data()

    result = data_np / flat_np

    newimage = dp.make_image(result, image, {})
    return newimage


# SPCAM keywords that should be added to the primary HDU
prihdr_kwds = [
    'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND', 'BZERO',
    'BSCALE', 'BUNIT', 'BLANK', 'DATE-OBS', 'UT', 'UT-STR', 'UT-END',
    'HST', 'HST-STR', 'HST-END', 'LST', 'LST-STR', 'LST-END', 'MJD',
    'TIMESYS', 'MJD-STR', 'MJD-END', 'ZD-STR', 'ZD-END', 'SECZ-STR',
    'SECZ-END', 'AIRMASS', 'AZIMUTH', 'ALTITUDE', 'PROP-ID', 'OBSERVER',
    'EXP-ID', 'DATASET', 'OBS-MOD', 'OBS-ALOC', 'DATA-TYP', 'OBJECT',
    'RA', 'DEC', 'RA2000', 'DEC2000', 'OBSERVAT', 'TELESCOP', 'FOC-POS',
    'TELFOCUS', 'FOC-VAL', 'FILTER01', 'EXPTIME', 'INSTRUME', 'INS-VER',
    'WEATHER', 'SEEING', 'ADC-TYPE', 'ADC-STR', 'ADC-END', 'INR-STR',
    'INR-END', 'DOM-WND', 'OUT-WND', 'DOM-TMP', 'OUT-TMP', 'DOM-HUM',
    'OUT-HUM', 'DOM-PRS', 'OUT-PRS', 'EXP1TIME', 'COADD', 'M2-POS1',
    'M2-POS2', 'M2-POS3', 'M2-ANG1', 'M2-ANG2', 'M2-ANG3', 'AUTOGUID',
    'COMMENT', 'INST-PA', 'EQUINOX',
    ]

# SPCAM keywords that should be added to the image HDUs
imghdr_kwds = [
    'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND', 'BZERO',
    'BSCALE', 'BUNIT', 'BLANK', 'FRAMEID', 'EXP-ID', 'DETECTOR', 'DET-ID',
    'DET-A01', 'DET-P101', 'DET-P201', 'DET-TMP', 'DET-TMED', 'DET-TMIN',
    'DET-TMAX', 'GAIN', 'EFP-MIN1', 'EFP-RNG1', 'EFP-MIN2', 'EFP-RNG2',
    'PRD-MIN1', 'PRD-RNG1', 'PRD-MIN2', 'PRD-RNG2', 'BIN-FCT1', 'BIN-FCT2',
    'DET-VER', 'S_UFNAME', 'S_FRMPOS', 'S_BCTAVE', 'S_BCTSD', 'S_AG-OBJ',
    'S_AG-RA', 'S_AG-DEC', 'S_AG-EQN', 'S_AG-X', 'S_AG-Y', 'S_AG-R',
    'S_AG-TH', 'S_ETMED', 'S_ETMAX', 'S_ETMIN', 'S_XFLIP', 'S_YFLIP',
    'S_M2OFF1', 'S_M2OFF2', 'S_M2OFF3', 'S_DELTAZ', 'S_DELTAD', 'S_SENT',
    'S_GAIN1', 'S_GAIN2', 'S_GAIN3', 'S_GAIN4', 'S_OSMN11', 'S_OSMX11',
    'S_OSMN21', 'S_OSMX21', 'S_OSMN31', 'S_OSMX31', 'S_OSMN41', 'S_OSMX41',
    'S_OSMN12', 'S_OSMX12', 'S_OSMN22', 'S_OSMX22', 'S_OSMN32', 'S_OSMX32',
    'S_OSMN42', 'S_OSMX42', 'S_EFMN11', 'S_EFMX11', 'S_EFMN21', 'S_EFMX21',
    'S_EFMN31', 'S_EFMX31', 'S_EFMN41', 'S_EFMX41', 'S_EFMN12', 'S_EFMX12',
    'S_EFMN22', 'S_EFMX22', 'S_EFMN32', 'S_EFMX32', 'S_EFMN42', 'S_EFMX42',
    'EQUINOX',  'CRVAL1',  'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2',
    'LONGPOLE', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2', 'WCS-ORIG',
    'RADECSYS', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
    ]

def spcam_pack(dirpath, exp_id, outpath='.', compress=False):
    """
    Pack a group of 10 SPCAM FITS files (a single exposure) into one
    FITS file with one primary HDU with no data and 10 image HDUs.

    Parameters
    ----------
    dirpath : str
        specifies the directory containing the exposure files
    exp_id : str
        a frame id of the format "SUPANNNNNNNN"--any one of the frame ids
        making up the single exposure
    outpath : str (optional)
        specifies the output directory for the resulting fits file.
        Defaults to current directory.
    compress : bool (optional)
        if True, will try to Rice-compress the image HDUs.  Note that
        this slows down the process considerably.  Default: False
    """

    import astropy.io.fits as pyfits

    exp_id = exp_id[:11]
    glob_exp = os.path.join(dirpath, exp_id + '?.fits')
    filelist = glob.glob(glob_exp)
    filelist.sort()

    hdus = {}

    outfile = exp_id.replace('A', 'E') + '0.fits'

    fitsobj = pyfits.HDUList()

    i = 0
    for fitsfile in filelist:
        with pyfits.open(fitsfile, 'readonly') as in_f:
            if i == 0:
                # prepare primary HDU header
                hdu = pyfits.PrimaryHDU()
                prihdr = hdu.header
                for kwd, val in in_f[0].header.items():
                    if kwd.upper() in prihdr_kwds:
                        #prihdr[kwd] = (val, comment)
                        prihdr[kwd] = val

                fitsobj.append(hdu)

            # create each image HDU
            data = numpy.copy(in_f[0].data.astype(numpy.uint16))

            if not compress:
                hdu = pyfits.ImageHDU(data=data)
            else:
                hdu = pyfits.CompImageHDU(data=data,
                                          compression_type='RICE_1')
            
            for kwd, val in in_f[0].header.items():
                if kwd.upper() in imghdr_kwds:
                    #hdu.header[kwd] = (val, comment)
                    hdu.header[kwd] = val

            # add CHECKSUM and DATASUM keywords
            if not compress:
                hdu.add_checksum()

            det_id = int(hdu.header['DET-ID'])
            #fitsobj.append(hdu)
            hdus[det_id] = hdu
            i += 1

    # stack HDUs in detector ID order
    for i in xrange(10):
        fitsobj.append(hdus[i])
        
    # fix up to FITS standard as much as possible
    fitsobj.verify('silentfix')

    # save to a file
    outpath = os.path.join(outpath, outfile)

    fitsobj.writeto(outpath)

#END
