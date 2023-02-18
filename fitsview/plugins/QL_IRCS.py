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
to and press the "Set Memo" button.  Multiple selection follows the usual
rules about holding down CTRL and/or SHIFT keys.

**Displaying an image**

Double-click on a log entry.

"""
import os
from collections import OrderedDict

from ginga import GingaPlugin, AstroImage
from ginga.misc import Bunch
from ginga.gw import Widgets

from g2base.astro.frame import Frame


__all__ = ['QL_IRCS']

import ObsLog

class QL_IRCS(ObsLog.ObsLog):

    def __init__(self, fv):
        super(QL_IRCS, self).__init__(fv)

        self.chname = 'IRCS'
        self.chnames = ['IRCS_Norm_Cam', 'IRCS_Norm_Spg']
        self.file_prefixes = ['IRCA']

        # columns to be shown in the table
        column_info = [dict(col_title="Array", fits_kwd='DET-ID'),
                       dict(col_title="Obs Mod", fits_kwd='OBS-MOD'),
                       dict(col_title="Datatype", fits_kwd='DATA-TYP'),
                       dict(col_title="FrameID", fits_kwd='FRAMEID',
                            col_width=160),
                       dict(col_title="Object", fits_kwd='OBJECT'),
                       #dict(col_title="UT", fits_kwd='UT'),
                       dict(col_title="HST", fits_kwd='HST'),
                       #dict(col_title="PropId", fits_kwd='PROP-ID'),
                       dict(col_title="Exp Time", fits_kwd='EXP1TIME'),
                       dict(col_title="Ndr", fits_kwd='NDR'),
                       dict(col_title="CoAdd", fits_kwd='COADD'),
                       dict(col_title="Air Mass", fits_kwd='AIRMASS'),
                       #dict(col_title="Pos Ang", fits_kwd='INST-PA'),
                       #dict(col_title="Ins Rot", fits_kwd='INSROT'),
                       #dict(col_title="Foc Val", fits_kwd='FOC-VAL'),
                       #dict(col_title="Filter01", fits_kwd='FILTER01'),
                       #dict(col_title="Filter02", fits_kwd='FILTER02'),
                       #dict(col_title="Filter03", fits_kwd='FILTER03'),
                       #dict(col_title="RA", fits_kwd='RA'),
                       #dict(col_title="DEC", fits_kwd='DEC'),
                       #dict(col_title="EQUINOX", fits_kwd='EQUINOX'),
                       dict(col_title="IMR STAT", fits_kwd='D_IMR'),
                       dict(col_title="PA", fits_kwd='D_IMRPAD'),
                       dict(col_title="IMR Mode", fits_kwd='D_IMRMOD'),
                       dict(col_title="CW1", fits_kwd='I_MCW1NM'),
                       dict(col_title="CW2", fits_kwd='I_MCW2NM'),
                       dict(col_title="CW3", fits_kwd='I_MCW3NM'),
                       dict(col_title="Cam Res", fits_kwd='I_CAMRES'),
                       dict(col_title="Cam Focus", fits_kwd='I_MFOCMC'),
                       dict(col_title="SLW", fits_kwd='I_SLWNM'),
                       dict(col_title="SPW", fits_kwd='I_SPWNM'),
                       dict(col_title="ECH", fits_kwd='I_MECHAS'),
                       dict(col_title="XDS", fits_kwd='I_MXDSAS'),
                       dict(col_title="Loop (AO)", fits_kwd='D_LOOP'),
                       dict(col_title="AO Mode", fits_kwd='D_MODE'),
                       dict(col_title="VM", fits_kwd='D_VMVOLT'),
                       dict(col_title="DM", fits_kwd='D_DMGAIN'),
                       dict(col_title="HTTG", fits_kwd='D_WTTG'),
                       dict(col_title="LTTG", fits_kwd='D_LTTG'),
                       dict(col_title="Memo", fits_kwd='G_MEMO'),
                       ]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QL_IRCS')
        self.settings.set(sortable=True,
                          color_alternate_rows=True,
                          column_info=column_info,
                          cache_normalized_images=True)
        self.settings.load(onError='silent')


    def build_gui(self, container):
        super(QL_IRCS, self).build_gui(container)

        self.w.obslog_dir.set_text("{}/Procedure/IRCS".format(os.environ['HOME']))

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
        if chname != 'IRCS' or not self.gui_up:
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
        fr = Frame(frameid)

        # normalized image prefix
        fr.frametype = 'N'
        newname = fr.frameid

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
        coadds = header.get('COADD', 1)
        ndr = header.get('NDR', 1)
        divisor = coadds * ndr

        data_np = image.get_data()
        data_np = data_np / divisor

        # create a new image
        new_image = AstroImage.AstroImage(data_np=data_np, logger=self.logger)
        new_image.set(name=newname)
        new_image.update_keywords(header)
        new_image.update_keywords(dict(COADD=1, NDR=1,
                                       FRAMEID=newname))

        if self.settings.get('cache_normalized_images', True):
            # write out a cached copy so we can reload as necessary
            try:
                prefix = os.path.join(os.environ['GEN2COMMON'],
                                      'data_cache', 'IRCS')
            except KeyError:
                prefix = '/tmp'
            cached_path = os.path.join(prefix, newname + '.fits')
            new_image.set(path=cached_path)
            if not os.path.exists(cached_path):
                new_image.save_as_file(cached_path)

        return new_image

    def view_image(self, frameid, info):
        chname = self.chnames[0] if info['DET-ID'] == 'CAM' else self.chnames[1]
        channel = self.fv.get_current_channel()
        if channel.name != chname:
            channel = self.fv.get_channel_on_demand(chname)
            self.fv.change_channel(chname)

        # want to see the normalized image
        imname = 'IRCN' + frameid[4:]
        if imname in channel:
            channel.switch_name(imname)

        else:
            #<-- need to load the original image and reprocess it
            chname = 'IRCS'
            # TODO: record the absolute path to the file in the ObsLog
            filepath = os.path.join('/gen2', 'share', 'data', chname,
                                    frameid + '.fits')
            self.logger.info(f"attempting to load '{filepath}'...")
            self.fv.load_file(filepath, chname=chname)

    def __str__(self):
        return 'ql_ircs'
