#
# focas.py -- FOCAS data processing routines
#
# E. Jeschke
#
#
"""
Dear Eric-san,

Thank you for the email.
The coordinates of the regions depend on CCD chip and binning.
The following numbers are x-coordinates starting from x=1 (not 0)
and BIN-FCT1=1,2,4 from left to right, respectively.

Chip1 (DET-ID=1)
  ch1 : 9-520, 5-260, 3-130
  ov1 : 521-536, 261-276, 131-146
  ov2 : 537-552, 277-292, 147-162
  ch2 : 553-1064, 293-548, 163-290
  ch3 : 1081-1592, 557-812, 295-422
  ov3 : 1593-1608, 813-828, 423-438
  ov4 : 1610-1625, 830-845, 440-455
  ch4 : 1626-2137, 846-1101, 456-583

Chip2 (DET-ID=2)
  ch1 : 9-520, 5-260, 3-130
  ov1 : 521-536, 261-276, 131-146
  ov2 : 537-552, 277-292, 147-162
  ch2 : 553-1064, 293-548, 163-290
  ch3 : 1081-1592, 557-812, 295-422
  ov3 : 1593-1608, 813-828, 423-438
  ov4 : 1609-1624, 829-844, 439-454
  ch4 : 1625-2136, 845-1100, 455-582

Region ov1,2,3,4 are used to calculate bias levels of ch1,2,3,4 and
ch1,2,3,4 are combined to create the overscan-region-removed image.
I prefer not to use the first and last pixels of an overscan region.
For example, in the case of BIN-FCT1=1, x=522-535 should be
used to subtract bias level from x=9-520.

For the WCS, it is correct.
CRPIX1,2 should be divided by BIN-FCT1,2.

I just realized that I haven't thought about partial readout mode.
It seems to require some investigation and modification of FOCAS
control program and let me get back to you later about this issue.
Anyway, dividing by BIN-FCT1,2 should work for full readout mode.

Thank you,
Takashi Hattori
"""
import os
import re, glob
import numpy

from ginga import AstroImage
from ginga.misc import Bunch, log
from ginga.util import dp, mosaic, io_fits

from g2base.astro.frame import Frame


class FocasDR(object):

    def __init__(self, logger=None):
        super(FocasDR, self).__init__()

        if logger is None:
            logger = log.get_logger(level=20, log_stderr=True)
        self.logger = logger

        self.num_ccds = 2
        self.num_frames = 2
        self.inscode = 'FCS'
        self.fov = 0.15

    def get_regions(self, image):
        """Extract the keywords defining the overscan and effective pixel
        regions in a FOCAS image.  The data is returned in a dictionary of
        bunches.  The keys of the dictionary are the channel numbers, plus
        'image'.
        """
        det_id = int(image.get_keyword("DET-ID"))   # {1, 2}
        binning = int(image.get_keyword("BIN-FCT1"))   # {1, 2, 4}

        # FOCAS images need to have CRPIX corrected for binning
        crpix1, crpix2 = [float(x)
                          for x in image.get_keywords_list("CRPIX1",
                                                           "CRPIX2"))]
        chip_d = focas_ccd['chip%d' % (det_id)]

        wd, ht = image.get_size()
        d = {}
        xcut = 0
        newwd = 0
        l = []
        for channel in (1, 2, 3, 4):
            efminx, efmaxx = chip_d['ch%d' % (channel)]['bin%d' % (binning)]
            # subtract one because these are 1 based
            efminx, efmaxx = efminx - 1, efmaxx - 1
            efminy, efmaxy = 0, ht - 1
            osminx, osmaxx = chip_d['ov%d' % (channel)]['bin%d' % (binning)]
            # subtract one because these are 1 based
            osminx, osmaxx = osminx - 1, osmaxx - 1
            osminy, osmaxy = 0, ht - 1
            xcut += osmaxx - osminx + 1
            newwd += efmaxx + 1 - efminx

            gain = float(image.get_keyword("GAIN"))
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
        print("width reduction=%d xcut=%d" % (wd - newwd, xcut))
        d['image'] = Bunch.Bunch(xcut=xcut, ycut=ycut,
                                 newwd=newwd, newht=newht)

        image.set_keyword('CRPIX1', crpix1 / binning)
        image.set_keyword('CRPIX2', crpix2 / binning)

        return d


    def subtract_overscan_np(self, data_np, d, header=None):
        """Subtract the median bias calculated from the overscan regions
        from a FOCAS image data array.  The resulting image is trimmed to
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
        if header is not None:
            header['NAXIS1'] = newwd
            header['NAXIS2'] = newht

        # original image size
        ht, wd = data_np.shape[:2]

        for channel in (1, 2, 3, 4):
            #print "processing channel %d" % (channel)
            ch = d[channel]

            # Hattori-san prefers to avoid the first and last pixel of
            # the overscan region for computing bias
            osminx, osmaxx = ch.osminx + 1, ch.osmaxx - 1

            # get median of each row in overscan area for this channel
            ovsc_median = numpy.median(data_np[ch.efminy:ch.efmaxy+1,
                                               osminx:osmaxx+1], axis=1)
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
        if flat_norm is not None:
            flat = flat / flat_norm

        img_flat = dp.make_image(flat, image, header)
        return img_flat


    def make_flat_tiles(self, datadir, explist, output_pfx='flat',
                        output_dir=None):

        # Get the median values for each CCD image
        flats = []
        for i in range(self.num_frames):
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

            if output_dir is None:
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
        for i in range(num_exp):
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


# NOTE: all indexes are 1-based!
focas_ccd = dict(
    chip1 = dict(
        ch1 = dict(bin1=(9, 520), bin2=(5, 260), bin4=(3, 130)),
        ov1 = dict(bin1=(521, 536), bin2=(261, 276), bin4=(131, 146)),
        ch2 = dict(bin1=(553, 1064), bin2=(293, 548), bin4=(163, 290)),
        ov2 = dict(bin1=(537, 552), bin2=(277, 292), bin4=(147, 162)),
        ch3 = dict(bin1=(1081, 1592), bin2=(557, 812), bin4=(295, 422)),
        ov3 = dict(bin1=(1593, 1608), bin2=(813, 828), bin4=(423, 438)),
        ch4 = dict(bin1=(1626, 2137), bin2=(846, 1101), bin4=(456, 583)),
        ov4 = dict(bin1=(1610, 1625), bin2=(830, 845), bin4=(440, 455)),
    ),
    chip2 = dict(
        ch1 = dict(bin1=(9, 520), bin2=(5, 260), bin4=(3, 130)),
        ov1 = dict(bin1=(521, 536), bin2=(261, 276), bin4=(131, 146)),
        ch2 = dict(bin1=(553, 1064), bin2=(293, 548), bin4=(163, 290)),
        ov2 = dict(bin1=(537, 552), bin2=(277, 292), bin4=(147, 162)),
        ch3 = dict(bin1=(1081, 1592), bin2=(557, 812), bin4=(295, 422)),
        ov3 = dict(bin1=(1593, 1608), bin2=(813, 828), bin4=(423, 438)),
        ch4 = dict(bin1=(1625, 2136), bin2=(845, 1100), bin4=(455, 582)),
        ov4 = dict(bin1=(1609, 1624), bin2=(829, 844), bin4=(439, 454)),
    ),
    )

#END
