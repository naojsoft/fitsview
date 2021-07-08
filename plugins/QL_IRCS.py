# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
"""IRCS Ginga QuickLook plugin.

**Plugin Type: Global**

``QL_IRCS`` is a global plugin. Only one instance can be opened.

**Usage**

**Saving the log to a file**

Put in values for the Observation Log folder and filename.  The format
of the file saved will depend on the file extension of the filename;
use the type selector combobox to pick the right extension:

* txt: whitespace separated, ascii commented header
* fits: binary table in a FITS file
* xlsx: MS Excel file format

The file is rewritten out every time a new entry is added to the log

**Adding a memo to one or more log entries**

Write a memo in the memo box.  Select one or more frames to add the memo
to and press the "Add Memo" button.  Multiple selection follows the usual
rules about holding down CTRL and/or SHIFT keys.

**Displaying an image**

Double-click on a log entry.

"""
import os
from collections import OrderedDict

from ginga import GingaPlugin, AstroImage
from ginga.misc import Bunch
from ginga.gw import Widgets

__all__ = ['QL_IRCS']

import ObsLog

class QL_IRCS(ObsLog.ObsLog):

    def __init__(self, fv):
        super(QL_IRCS, self).__init__(fv)

        self.chname = 'IRCS'
        self.chnames = ['IRCS_Norm_Cam', 'IRCS_Norm_Spg']
        self.file_prefixes = ['IRCA']
        #self.ql_tagname = 'AQUISITION'

        # columns to be shown in the table
        columns = [("Array", 'DET-ID'),
                   ("Obs Mod", 'OBS-MOD'),
                   ("Datatype", 'DATA-TYP'),
                   ("FrameID", 'FRAMEID'),
                   ("Object", 'OBJECT'),
                   #("UT", 'UT'),
                   ("HST", 'HST'),
                   #("PropId", 'PROP-ID'),
                   ("Exp Time", 'EXPTIME'),
                   ("Ndr", 'NDR'),
                   ("CoAdds", 'COADDS'),
                   ("Air Mass", 'AIRMASS'),
                   #("Pos Ang", 'INST-PA'),
                   #("Ins Rot", 'INSROT'),
                   #("Foc Val", 'FOC-VAL'),
                   #("Filter01", 'FILTER01'),
                   #("Filter02", 'FILTER02'),
                   #("Filter03", 'FILTER03'),
                   #("RA", 'RA'),
                   #("DEC", 'DEC'),
                   #("EQUINOX", 'EQUINOX'),
                   ("IMR STAT", 'D_IMR'),
                   ("PA", 'D_IMRPAD'),
                   ("IMR Mode", 'D_IMRMOD'),
                   ("CW1", 'I_MCW1NM'),
                   ("CW2", 'I_MCW2NM'),
                   ("CW3", 'I_MCW3NM'),
                   ("Cam Res", 'I_CAMRES'),
                   ("Cam Focus", 'I_MFOCMC'),
                   ("SLW", 'I_SLWNM'),
                   ("SPW", 'I_SPWNM'),
                   ("ECH", 'I_MECHAS'),
                   ("XDS", 'I_MXDSAS'),
                   ("Loop (AO)", 'D_LOOP'),
                   ("AO Mode", 'D_MODE'),
                   ("VM", 'D_VMVOLT'),
                   ("DM", 'D_DMGAIN'),
                   ("HTTG", 'D_WTTG'),
                   ("LTTG", 'D_LTTG'),
                   ("Memo", 'G_MEMO'),
                  ]

        self.settings.set(sortable=True,
                          color_alternate_rows=True,
                          report_columns=columns,
                          cache_normalized_images=True)

    def replace_kwds(self, header):
        d = super(QL_IRCS, self).replace_kwds(header)

        d['DET-ID'] = 'CAM' if str(d.get('DET-ID', 1)).strip() == '1' else 'SPG'

        if d.get('D_IMR', '').upper().strip() != 'TRACK':
            d['D_IMRPAD'] = "---"
            d['D_IMRMOD'] = "---"

        if d.get('D_LOOP', '').upper().strip() == 'ON':
            if d.get('D_MODE', '').upper().strip() != 'LGS':
                d['D_WTTG'] = "---"
                d['D_LTTG'] = "---"

        else:
            d['D_MODE'] = "---"
            d['D_VMVOLT'] = "---"
            d['D_DMGAIN'] = "---"
            d['D_WTTG'] = "---"
            d['D_LTTG'] = "---"

        return d

    def process_image(self, chname, header, image):
        if chname != 'IRCS':
            return

        imname = image.get('name', None)
        if imname is None:
            return

        # skip pre-baked "quicklook" (Q) frames
        if imname.startswith('IRCQ'):
            return

        header = image.get_header()

        frameid = header.get('FRAMEID', None)
        if frameid is None:
            return

        frameid = frameid.strip()

        # normalized image prefix
        newname = 'IRCN' + frameid[4:]

        try:
            det_id = int(header.get('DET-ID', 1)) - 1
        except Exception as e:
            self.logger.error("Error getting DET-ID: {}".format(e))
            det_id = 0
        chname = self.chnames[det_id]

        channel = self.fv.get_channel_on_demand(chname)

        # check if the image has already been processed
        if newname in channel:
            return

        # otherwise make a normalized image and add it to our channel
        new_image = self.make_normalized_image(newname, image)

        self.fv.gui_do(channel.add_image, new_image)

    def make_normalized_image(self, newname, image):
        header = image.get_header()

        # normalize the data
        coadds = header.get('COADDS', 1)
        ndr = header.get('NDR', 1)
        divisor = coadds * ndr

        data_np = image.get_data()
        data_np = data_np / divisor

        # create a new image
        new_image = AstroImage.AstroImage(data_np=data_np, logger=self.logger)
        new_image.set(name=newname)
        new_image.update_keywords(header)
        new_image.update_keywords(dict(COADD=1, COADDS=1, NDR=1,
                                       FRAMEID=newname))

        if self.settings.get('cache_normalized_images', True):
            # write out a cached copy so we can reload as necessary
            try:
                prefix = os.path.join(os.environ['GEN2COMMON'], 'data_cache')
            except KeyError:
                prefix = '/tmp'
            cached_path = os.path.join(prefix, newname + '.fits')
            new_image.set(path=cached_path)
            if not os.path.exists(cached_path):
                new_image.save_as_file(cached_path)

        return new_image

    def view_image(self, imname, info):
        chname = self.chnames[0] if info['DET-ID'] == 'CAM' else self.chnames[1]
        channel = self.fv.get_current_channel()
        if channel.name != chname:
            channel = self.fv.get_channel(chname)
            self.fv.change_channel(chname)

        # want to see the normalized image
        imname = 'IRCN' + imname[4:]
        channel.switch_name(imname)

    def __str__(self):
        return 'ql_ircs'
