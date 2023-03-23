#
# fitsUtils.py -- A utility file with methods to manipulate FITS files
# Works in conjunction with MESOffset ginga plugin for MOS Acquisition
#
#   2016/Aug/01  Justin Kunimune   Initial version
#
#   2017/Sep/15  Dr. Chi-Hung Yan (chyan@naoj.org, chyan@asiaa.sinica.edu.tw)
#
#      Adding functions to remove the denpendicies of IRAF.
#
#   2017/Sep/12  Dr. Chi-Hung Yan (chyan@naoj.org, chyan@asiaa.sinica.edu.tw)
#
#      Adding debug mode for stopping produce by-products.
#


# standard imports
import os

# third-party imports
from astropy.io import fits
import astropy.utils.introspection
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import shift

# The astropy.io.fits.writeto API changed from a "clobber" argument to
# an "overwrite" argument in version 1.3. Determine which version we
# are running and set the keyword argument appropriately.
if astropy.utils.introspection.minversion(astropy, '1.3'):
    write_to_kwargs = {'overwrite': True}
else:
    write_to_kwargs = {'clobber': True}

# constants
# DIR_MCSRED = '../../MCSRED2/'
DIR_MCSRED = os.path.join(os.environ['HOME'], 'Procedure/MOIRCS/MCSRED2/')
NO_SUCH_FILE_ERR = ("No such file or directory: {}\nPlease check your frame " +
                    "numbers and image directory, or run Ginga from a " +
                    "different directory.")
WRONG_CHIP_ERR = ("{} should be data from chip {}, but is from chip {}. Try " +
                  "a different frame number.")
LOW_ELEV_WARN = (u"{}MCSA{:08d}.fits has low elevation of {:.1f}\u00B0; the " +
                 "mosaicing database may not be applicable here.")
USER_INTERRUPT_ERR = ("This process was terminated. Please press 'Return to " +
                      "Menu' to start it over.")

SAVE_INTERMEDIATE_FILES = True


# def iraf_login_cl():
#     # To be able to import the pyraf modules, there has to be a
#     # login.cl file in either $HOME or $HOME/.iraf. If there isn't,
#     # create $HOME/login.cl and put tell IRAF to load the "image"
#     # operators, which includes geotran and imcombine.
#     login_cl_filenames = [os.path.join(os.environ['HOME'], 'login.cl'),
#                           os.path.join(os.environ['HOME'], '.iraf',
#                           'login.cl')]
#     file_exists = False
#     for filename in login_cl_filenames:
#         if os.path.exists(filename):
#             file_exists = True
#     if not file_exists:
#         with open(login_cl_filenames[0], 'w') as outputfile:
#             outputfile.write('images  # general image operators\n')


# iraf_login_cl()
# from pyraf.iraf import geotran
# from pyraf.iraf import imcombine
# from pyraf.iraf import rotate


class DistortionCoeffCh1:
    a = -8.635143
    b = 1.009153
    c = 2.394803E-4
    d = -10.43176
    e = 2.192026E-4
    f = 1.00965

    xc00 = -16.37997
    xc10 = 0.04216358
    xc20 = -3.961311E-5
    xc30 = 1.453919E-8
    xc01 = 0.02687643
    xc11 = -2.919222E-5
    xc21 = -6.896180E-11
    xc02 = -1.334916E-5
    xc12 = 1.432681E-8
    xc03 = 6.087694E-11

    yc00 = -16.98225
    yc10 = 0.02702253
    yc20 = -1.542926E-5
    yc30 = 2.064173E-10
    yc01 = 0.04707864
    yc11 = -2.632263E-5
    yc21 = 1.465591E-8
    yc02 = -4.374318E-5
    yc12 = -2.114975E-12
    yc03 = 1.430626E-8


class DistortionCoeffCh2:
    a = -5.023962
    b = 1.004285
    c = 1.277659E-4
    d = -5.48291
    e = 7.625483E-5
    f = 1.004831

    xc00 = -27.69665
    xc10 = 0.06082649
    xc20 = -4.797236E-5
    xc30 = 1.420412E-8
    xc01 = 0.03267656
    xc11 = -2.884639E-5
    xc21 = -9.848598E-11
    xc02 = -1.614414E-5
    xc12 = 1.410305E-8
    xc03 = 9.066495E-11

    yc00 = -24.74531
    yc10 = 0.03398627
    yc20 = -1.568090E-5
    yc30 = 2.707538E-10
    yc01 = 0.05459096
    yc11 = -3.218881E-5
    yc21 = 1.433751E-8
    yc02 = -4.413517E-5
    yc12 = -4.351142E-11
    yc03 = 1.435618E-8


class MosaicParameterCh2:
    a = -1602.189
    b = 0.9993759
    c = 0.009615554
    d = 40.30523
    e = -0.009611026
    f = 0.9998468


def nothing(*args, **kwargs):
    """A placeholder function for log"""
    pass


def auto_process_fits(mode, n1, n2, c, i, w, f, t,
                      log=nothing, next_step=None):
    """
    Use mode to choose a fitsUtils method and call it with the appropriate
    arguments
    @param mode:
        A string which will determine the method we call - either 'star',
        'mask', or 'starhole'
    """
    try:
        if mode == 'mask':
            process_mask_fits(n1, c, i, w, f, t, log, next_step=next_step)
        else:
            process_star_fits(n1, n2, c, i, w, f, t, log, next_step=next_step)
        if t.is_set():
            raise RuntimeError(USER_INTERRUPT_ERR)
    except Exception as e:
        log("{}: {}".format(type(e).__name__, e), level='e')


def process_star_fits(star_num, back_num, c_file, img_dir, work_dir,
                      output_filename,
                      terminate, log=nothing, next_step=None):
    """
    Process the raw star and background images by subtracting the background
    from the star images, adding a gaussian filter to the result, and mosaicing
    it all together
    @param star_num:
        The integer in the star image chip1 filename
    @param back_num:
        The integer in the background image chip1 filename
    @param c_file:
        The location of the cfg file that contains parameters for make_mosaic
    @param img_dir:
        The string prefix to all raw image filenames
    @param work_dir:
        The string prefix to all intermediate processing filenames
    @param output_filename:
        The filename of the final FITS image
    @param terminate:
        The threading.Event object that flags True if the user tries to
        terminate this thread
    @param log:
        A function which should take one argument, and will be called to report
        information
    @param next_step:
        The function to be called at the end of this process
    @raises IOError:
        If it cannot find the FITS files in the specified directory
    @raises ValueError:
        If the FITS files have the wrong chip values
    """
    log("Processing star frames...")

    # open the star FITS and check header info
    star_chip = []
    for i in (0, 1):
        star_chip.append(open_fits("{}MCSA{:08d}.fits".format(
            img_dir, star_num + i), i + 1))
        if star_chip[i].header['ALTITUDE'] < 45.0:
            log(LOW_ELEV_WARN.format(img_dir, star_num + i,
                                     star_chip[i].header['ALTITUDE']),
                level='warning')

    # subtract the background frames from the star frames
    log("Subtracting images...")
    if back_num != 0:
        back_chip = []
        for i in (0, 1):
            back_chip.append(open_fits("{}MCSA{:08d}.fits".format(
                img_dir, back_num + i), i + 1))

        dif_data = [star_chip[i].data - back_chip[i].data for i in (0, 1)]

    else:
        dif_data = [star_chip[i].data for i in (0, 1)]

    # mosaic the chips together
    if terminate.is_set():
        return
    mosaic_data = makeMosaic(
        'star',
        star_num,
        dif_data,
        c_file,
        work_dir,
        terminate,
        log=log)
    # mosaic_data = make_mosaic('star', star_num, dif_data, \
    #  c_file, work_dir, terminate, log=log)
    if terminate.is_set():
        return

    # apply gaussian blur
    log("Blurring...")
    mosaic_data = gaussian_filter(mosaic_data, 1.0)

    # write to file and go to next_step
    fits.writeto(output_filename, mosaic_data, header=star_chip[0].header,
                 **write_to_kwargs)
    if next_step is not None:
        next_step()


def process_mask_fits(mask_num, c_file, img_dir, work_dir, output_filename,
                      terminate, log=nothing, next_step=None):
    """
    Process the raw mask frames by changing their data type and mosaicing them
    together
    @param mask_num:
        The number in the filename of the mask chip1 FITS image
    @param c_file:
        The location of the cfg file that controls make_mosaic
    @param img_dir:
        The prefix for all raw image filenames
    @param work_dir:
        The string prefix to all intermediate processing filenames
    @param output_filename:
        The filename of the output FITS image
    @param terminate:
        The threading.Event that will flag True if the user tries to kill this
        process
    @param log:
        The function that will be called whenever something interesting happens
    @param next_step:
        The function to be called at the end of this process
    @raises IOError:
        If it cannot find the FITS images
    """
    log("Processing mask frames...")

    # load the files
    mask_chip = []
    for chip in (0, 1):
        mask_chip.append(open_fits("{}MCSA{:08d}.fits".format(
            img_dir, mask_num + chip), chip + 1))

    # mosaic the reformatted results to a file
    if terminate.is_set():
        return
    mosaic_data = makeMosaic('mask', mask_num,
                             [hdu.data for hdu in mask_chip], c_file,
                             work_dir, terminate, log=log)
    # mosaic_data = make_mosaic('mask', mask_num,
    #                          [hdu.data for hdu in mask_chip], c_file,
    #                          work_dir, terminate, log=log)
    if terminate.is_set():
        return

    # finish up by writing to file and moving on
    fits.writeto(output_filename, mosaic_data, header=mask_chip[0].header,
                 **write_to_kwargs)
    if next_step is not None:
        next_step()


def transformLocation(uncorrectXY, dc):

    x, y = uncorrectXY

    xfit1 = dc.a + dc.b * x + dc.c * y
    yfit1 = dc.d + dc.e * x + dc.f * y

    xfit2 = dc.xc00 + \
        (dc.xc10 * x) + \
        (dc.xc01 * y) + \
        (dc.xc20 * x**2) + \
        (dc.xc30 * x**3) + \
        (dc.xc11 * x * y) + \
        (dc.xc21 * x**2 * y) + \
        (dc.xc02 * y**2) + \
        (dc.xc12 * x * y**2) + \
        (dc.xc03 * y**3)

    yfit2 = dc.yc00 + \
        (dc.yc10 * x) + \
        (dc.yc20 * x**2) + \
        (dc.yc30 * x**3) + \
        (dc.yc01 * y) + \
        (dc.yc11 * x * y) + \
        (dc.yc21 * x**2 * y) + \
        (dc.yc02 * y**2) + \
        (dc.yc12 * x * y**2) + \
        (dc.yc03 * y**3)

    xfit = xfit1 + xfit2
    yfit = yfit1 + yfit2

    xfit[np.where(xfit < 0)] = 0
    yfit[np.where(yfit < 0)] = 0
    xfit[np.where(xfit > 2047)] = 2047
    yfit[np.where(yfit > 2047)] = 2047

    indx = np.uint(np.round(xfit))
    indy = np.uint(np.round(yfit))

    return indx, indy


def transformImage(base_name, n, input_arr, dc):
    """
    Correct the input array for distortion using the given dbs and gmp
    @param input_arr:
        The input numpy array
    @param dbs_filename:
        The filename of the IRAF 'database file'. Not to be confused with an
        SQLBase .dbs file, or a .db database file.
    @param gmp_filename:
        The filename inside the filename that tells IRAF which part of the dbs
        file to look at. Because using integers was too easy, so why not just
        use filenames with an extention that doesn't exist to index through a
        file. Rather, the .gmp file extension does exist, and serves a log of
        purposes, but none of them have anything to do with IRAF or image
        transformations. WHYYYYYY?
    @returns:
        The corrected numpy array
    """
    input_filename = base_name + '_geotran_input%s.fits' % n
    output_filename = base_name + '_geotran_output%s.fits' % n

    if os.path.exists(output_filename):
        os.remove(output_filename)

    if SAVE_INTERMEDIATE_FILES is True:
        fits.PrimaryHDU(data=input_arr).writeto(input_filename, **write_to_kwargs)

    # remapImage=detrendImage
    img = input_arr

    y, x = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    indx, indy = transformLocation([x, y], dc)

    remap = img[indy, indx]

    if SAVE_INTERMEDIATE_FILES is True:
        fits.PrimaryHDU(data=remap).writeto(output_filename, **write_to_kwargs)

    # output = fits.open(output_filename)[0].data
    #  if not SAVE_INTERMEDIATE_FILES:
    #    os.remove(input_filename)
    #    os.remove(output_filename)

    return remap


def makeMosaic(im_type, frnum, input_data, c_file,
               work_dir, terminate, log=nothing):
    """
    Correct the images for distortion, and then combine the two FITS images by
    rotating and stacking them vertically. Also do something to the header
    @param input_data:
        A sequence of two numpy 2D arrays to mosaic together
    @param c_file:
        The location of the configuration .cfg file that manages distortion-
        correction
    @param work_dir:
        The string prefix to all intermediate processing filenames
    @param terminate:
        The threading.Event object that will tell us when/if to terminate
    @param log:
        A function that takes a single string argument and records it somehow
    @returns:
        A mosaiced numpy array comprising data from the two input_data arrays
    """
    # read MSCRED c_file
    cfg = open(c_file, 'r')
    config = []
    line = cfg.readline()
    while line != '':
        if line[0] != '#':
            config.append(line.split()[-1].replace('dir_mcsred$', DIR_MCSRED))
        line = cfg.readline()
    cfg.close()

    mosaic_data = [None, None]
    if terminate.is_set():
        return

    basename1 = os.path.join(work_dir, im_type + '_MCSA{:08d}'.format(frnum))
    basename2 = os.path.join(
        work_dir,
        im_type +
        '_MCSA{:08d}'.format(
            frnum +
            1))

    # correct for distortion and apply mask
    log("fitsUtil debug mode = %s" % SAVE_INTERMEDIATE_FILES)

    log("Correcting for distortion using Python...")
    dcc1 = DistortionCoeffCh1()
    mosaic_data[0] = transformImage(basename1, 1, input_data[0], dcc1)
    if terminate.is_set():
        return

    dcc2 = DistortionCoeffCh2()
    mosaic_data[1] = transformImage(basename2, 1, input_data[1], dcc2)
    if terminate.is_set():
        return

    log("Correcting for more distortion using Python...")
    ''' Shift the channel 1 image '''
    temp1 = np.zeros((2048, 3636))
    temp1[:, 0:2048] = mosaic_data[0]
    mosaic_data[0] = temp1
    if terminate.is_set():
        return

    '''Operate on channel2'''
    temp2 = np.zeros((2048, 3636))
    temp2[:, 0:2048] = mosaic_data[1]

    y, x = np.mgrid[0:temp2.shape[0], 0:temp2.shape[1]]

    '''Loading mosaic parameter'''
    mc2 = MosaicParameterCh2()
    xfit = mc2.a + mc2.b * x + mc2.c * y
    yfit = mc2.d + mc2.e * x + mc2.f * y

    ''' Arranging indexes in correct range'''
    xfit[np.where(xfit < 0)] = 0
    yfit[np.where(yfit < 0)] = 0
    xfit[np.where(xfit > 2047)] = 2047
    yfit[np.where(yfit > 2047)] = 2047

    ''' Round the float numbers to integer'''
    indx = np.uint(np.round(xfit))
    indy = np.uint(np.round(yfit))

    remap = temp2[indy, indx]
    mosaic_data[1] = remap
    if terminate.is_set():
        return
    # combine and rotate the images
    log("Combining the chips and applying bad pixel mask using Python...")
    combined_image = combineMask(basename1, basename2, mosaic_data[0],
                                 mosaic_data[1], config[12], config[13])
    if terminate.is_set():
        return

    log("Rotating the image with Python...")
    input_filename = basename1 + '_rotate_input.fits'
    output_filename = basename1 + '_rotate_output.fits'
    if os.path.exists(input_filename):
        os.remove(input_filename)

    if SAVE_INTERMEDIATE_FILES is True:
        fits.PrimaryHDU(
            data=combined_image).writeto(
            input_filename,
            **write_to_kwargs)

    mosaic_arr = shift(
        np.rot90(
            combined_image, k=3), [
            33.5, 0], order=0)[
                67:, :]

    # log("SAVE_INTERMEDIATE_FILES =... %r " % SAVE_INTERMEDIATE_FILES)
    if SAVE_INTERMEDIATE_FILES is True:
        fits.PrimaryHDU(data=mosaic_arr).writeto(output_filename,
                                                 **write_to_kwargs)

    # mosaic_arr = fits.open(output_filename)[0].data
    # if not SAVE_INTERMEDIATE_FILES:
    #    os.remove(input_filename)
    #    os.remove(output_filename)
    if terminate.is_set():
        return
    # XXX: stuff I haven't figured out how to do wiothout IRAF yet :XXX #

    return mosaic_arr


# def make_mosaic(im_type, frnum, input_data, c_file,
#                 work_dir, terminate, log=nothing):
#     """
#     Correct the images for distortion, and then combine the two FITS
#     images by rotating and stacking them vertically. Also do something
#     to the header
#     @param input_data:
#         A sequence of two numpy 2D arrays to mosaic together
#     @param c_file:
#         The location of the configuration .cfg file that manages distortion-
#         correction
#     @param work_dir:
#         The string prefix to all intermediate processing filenames
#     @param terminate:
#         The threading.Event object that will tell us when/if to terminate
#     @param log:
#         A function that takes a single string argument and records
#         it somehow
#     @returns:
#         A mosaiced numpy array comprising data from the two input_data
#         arrays
#     """
#     # read MSCRED c_file
#     cfg = open(c_file, 'r')
#     config = []
#     line = cfg.readline()
#     while line != '':
#         if line[0] != '#':
#             config.append(line.split()[-1].replace('dir_mcsred$',
#                           DIR_MCSRED))
#         line = cfg.readline()
#     cfg.close()

#     mosaic_data = [None, None]
#     if terminate.is_set():
#         return

#     basename1 = os.path.join(work_dir, im_type +
#                          '_MCSA{:08d}'.format(frnum))
#     basename2 = os.path.join(
#         work_dir,
#         im_type +
#         '_MCSA{:08d}'.format(
#             frnum +
#             1))

#     # XXX: stuff I haven't figured out how to do wiothout IRAF yet :XXX #
#     # correct for distortion and apply mask
#     log("Correcting for distortion...")
#     mosaic_data[0] = transform(
#         basename1,
#         1,
#         input_data[0],
#         config[2],
#         config[3])
#     if terminate.is_set():
#         return
#     mosaic_data[1] = transform(
#         basename2,
#         1,
#         input_data[1],
#         config[4],
#         config[5])
#     if terminate.is_set():
#         return

#     log("Correcting for more distortion...")
#     mosaic_data[0] = transform(
#         basename1,
#         2,
#         mosaic_data[0],
#         config[8],
#         config[9])
#     if terminate.is_set():
#         return
#     mosaic_data[1] = transform(
#         basename2,
#         2,
#         mosaic_data[1],
#         config[10],
#         config[11])
#     if terminate.is_set():
#         return

#     # combine and rotate the images
#     log("Combining the chips and applying bad pixel mask...")

#     combined_image = \
#       combine_and_apply_mask(
#             basename1,
#             basename2,
#             mosaic_data[0],
#             mosaic_data[1],
#            config[12],
#            config[13]
#            )

#     if terminate.is_set():
#         return

#     log("Rotating the image...")
#     input_filename = basename1 + '_rotate_input.fits'
#     output_filename = basename1 + '_rotate_output.fits'
#     if os.path.exists(input_filename):
#         os.remove(input_filename)
#     fits.PrimaryHDU(data=combined_image).writeto(input_filename,
#                                     **write_to_kwargs)
#     rotate(input_filename, output_filename, 90.0, ncols=2048, nlines=3569)
#     mosaic_arr = fits.open(output_filename)[0].data

#     if not SAVE_INTERMEDIATE_FILES:
#         os.remove(input_filename)
#         os.remove(output_filename)
#     if terminate.is_set():
#         return
#     # XXX: stuff I haven't figured out how to do wiothout IRAF yet :XXX #

#     return mosaic_arr


def open_fits(filename, chipnum):
    """
    It's like astropy.fits.open, but with better error handling
    @param filename:
        The name of the FITS file to be opened ('***.fits')
    @param chipnum:
        The desired value of the DET-ID header card
    @returns:
        An astropy HDU object -- the first in the FITS file
    @raises IOError:
        if the file cannot be found
    @raises ValueError:
        if the DET-ID is not chipnum
    """
    try:
        hdu = fits.open(filename, 'readonly', memmap=False)[0]
    except IOError as e:
        if len(filename) >= 1 and filename[0] in ('/'):
            raise IOError(NO_SUCH_FILE_ERR.format(filename))
        else:
            raise IOError(
                NO_SUCH_FILE_ERR.format(
                    os.getcwd() + "/" + filename))
    if hdu.header['DET-ID'] != chipnum:
        raise ValueError(WRONG_CHIP_ERR.format(filename, chipnum,
                                               hdu.header['DET-ID']))
    return hdu


# XXX: Methods that still use IRAF :XXX #

# def transform(base_name, n, input_arr, dbs_filename, gmp_filename):
#     """
#     Correct the input array for distortion using the given dbs and gmp
#     @param input_arr:
#         The input numpy array
#     @param dbs_filename:
#         The filename of the IRAF 'database file'. Not to be confused with
#         an SQLBase .dbs file, or a .db database file.
#     @param gmp_filename:
#         The filename inside the filename that tells IRAF which part of
#         the dbs file to look at. Because using integers was too easy,
#         so why not just use filenames with an extention that doesn't
#         exist to index through a file. Rather, the .gmp file extension
#         does exist, and serves a log of purposes, but none of them have
#         anything to do with IRAF or image transformations. WHYYYYYY?
#     @returns:
#         The corrected numpy array
#     """
#     input_filename = base_name + '_geotran_input%s.fits' % n
#     output_filename = base_name + '_geotran_output%s.fits' % n
#     if os.path.exists(output_filename):
#         os.remove(output_filename)
#     fits.PrimaryHDU(data=input_arr).writeto(input_filename, **write_to_kwargs)
#     #geotran(input_filename, output_filename, dbs_filename, gmp_filename,
#     #        verbose='yes')
#     output = fits.open(output_filename)[0].data
#     if not SAVE_INTERMEDIATE_FILES:
#         os.remove(input_filename)
#         os.remove(output_filename)
#     return output


# def combine_and_apply_mask(basename1, basename2, input1,
#                     input2, filename_mask1, filename_mask2, mask_val=0):
#     """
#     Replace all masked pixels with zero in the input array
#     @param input1:
#         The input numpy array for chip 1
#     @param input2:
#         The input numpy array for chip 2
#     @param filename_mask1:
#         Filename of pixel mask for chip 1
#     @param filename_mask2:
#         Filename of pixel mask for chip 2
#     @param mask_val:
#         The value to put into all of the masked pixels
#     """

#     input_filename1 = basename1 + '_imcombine_input.fits'
#     input_filename2 = basename2 + '_imcombine_input.fits'
#     output_filename = basename1 + '_imcombine_output.fits'
#     if os.path.exists(output_filename):
#         os.remove(output_filename)

#     hdu = fits.PrimaryHDU(data=input1)
#     hdu.header['BPM'] = filename_mask1
#     hdu.writeto(input_filename1, **write_to_kwargs)

#     hdu = fits.PrimaryHDU(data=input2)
#     hdu.header['BPM'] = filename_mask2
#     hdu.writeto(input_filename2, **write_to_kwargs)

#     imcombine(','.join([input_filename1, input_filename2]), output_filename,
#          combine='average', reject='avsig',
#          masktype='goodvalue', maskvalue=mask_val)

#     output = fits.open(output_filename)[0].data

#     if not SAVE_INTERMEDIATE_FILES:
#         os.remove(input_filename1)
#         os.remove(input_filename2)
#         os.remove(output_filename)

#     return output


def combineMask(basename1, basename2, input1, input2,
                filename_mask1, filename_mask2, mask_val=0):
    """
    Replace all masked pixels with zero in the input array
    @param input1:
        The input numpy array for chip 1
    @param input2:
        The input numpy array for chip 2
    @param filename_mask1:
        Filename of pixel mask for chip 1
    @param filename_mask2:
        Filename of pixel mask for chip 2
    @param mask_val:
        The value to put into all of the masked pixels
    """

    input_filename1 = basename1 + '_imcombine_input.fits'
    input_filename2 = basename2 + '_imcombine_input.fits'
    output_filename = basename1 + '_imcombine_output.fits'
    if os.path.exists(output_filename):
        os.remove(output_filename)

    hdu = fits.PrimaryHDU(data=input1)
    badpixfit1 = fits.open(filename_mask1 + '.fits')[0]
    badpix1 = badpixfit1.data

    if SAVE_INTERMEDIATE_FILES is True:
        hdu.writeto(input_filename1, **write_to_kwargs)

    hdu = fits.PrimaryHDU(data=input2)
    badpixfit2 = fits.open(filename_mask2 + '.fits')[0]
    badpix2 = badpixfit2.data
    if SAVE_INTERMEDIATE_FILES is True:
        hdu.writeto(input_filename2, **write_to_kwargs)

    img = input1 * np.subtract(1, badpix1) + input2 * np.subtract(1, badpix2)

    avg = np.median(img)
    img = np.add(0.5 * avg, 0.5 * img)

    outhdu = fits.PrimaryHDU(data=img)
    if SAVE_INTERMEDIATE_FILES is True:
        # log("Writing commbined image...")
        outhdu.writeto(output_filename, **write_to_kwargs)
    # output = fits.open(output_filename)[0].data

#    if not SAVE_INTERMEDIATE_FILES:
#        os.remove(input_filename1)
#        os.remove(input_filename2)
#        os.remove(output_filename)

    # Close images after using them
    badpixfit1._close()
    badpixfit2._close()

    return img

# END
