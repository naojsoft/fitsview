# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
"""MOIRCS Ginga QuickLook plugin.

**Plugin Type: Global**

``QL_MOIRCS`` is a global plugin. Only one instance can be opened.

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

from astropy.io import fits

from ginga import GingaPlugin, AstroImage
from ginga.misc import Bunch
from ginga.gw import Widgets

from g2base.astro.frame import Frame


__all__ = ['QL_MOIRCS']

import ObsLog

class QL_MOIRCS(ObsLog.ObsLog):

    def __init__(self, fv):
        super().__init__(fv)

        self.chnames = ['MOIRCS_1', 'MOIRCS_2']
        self.file_prefixes = ['MCSA']

        # columns to be shown in the table
        column_info = [#dict(col_title="Array", fits_kwd='DET-ID'),
                       #dict(col_title="Obs Mod", fits_kwd='OBS-MOD'),
                       #dict(col_title="Datatype", fits_kwd='DATA-TYP'),
                       dict(col_title=self.sort_hdr, fits_kwd='FRAMEID'),
                       dict(col_title="Date(UT)", fits_kwd='DATE-OBS'),
                       dict(col_title="Time(HST)", fits_kwd='HST-STR'),
                       #dict(col_title="PropId", fits_kwd='PROP-ID'),
                       dict(col_title="Object", fits_kwd='OBJECT'),
                       dict(col_title="Exp Time", fits_kwd='EXPTIME'),
                       dict(col_title="Det-NSMP", fits_kwd='DET-NSMP'),
                       dict(col_title="K_DITCNT", fits_kwd='K_DITCNT'),
                       dict(col_title="K_EXPCNT", fits_kwd='K_EXPCNT'),
                       dict(col_title="Filter01", fits_kwd='FILTER01'),
                       dict(col_title="Filter02", fits_kwd='FILTER02'),
                       dict(col_title="Filter03", fits_kwd='FILTER03'),
                       dict(col_title="Air Mass", fits_kwd='AIRMASS'),
                       dict(col_title="SLIT", fits_kwd='SLIT'),
                       #dict(col_title="UT", fits_kwd='UT'),
                       #dict(col_title="Pos Ang", fits_kwd='INST-PA'),
                       #dict(col_title="Ins Rot", fits_kwd='INSROT'),
                       #dict(col_title="Foc Val", fits_kwd='FOC-VAL'),
                       #dict(col_title="RA", fits_kwd='RA'),
                       #dict(col_title="DEC", fits_kwd='DEC'),
                       #dict(col_title="EQUINOX", fits_kwd='EQUINOX'),
                       dict(col_title="Memo", fits_kwd='G_MEMO'),
                       ]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QL_MOIRCS')
        self.settings.set(sortable=True,
                          color_alternate_rows=True,
                          column_info=column_info,
                          cache_normalized_images=True)
        self.settings.load(onError='silent')

        self.col_info = self.settings.get('column_info', [])
        # this will set rpt_columns and col_widths
        self.process_columns(self.col_info)

    def build_gui(self, container):
        super().build_gui(container)

        self.w.obslog_dir.set_text("{}/Procedure/MOIRCS".format(os.environ['HOME']))

        self.w.auto_save.set_state(True)

    def replace_kwds(self, header):
        d = super().replace_kwds(header)
        d['SLIT'] = d.get('SLIT', '---').strip()
        return d

    def process_image(self, chname, header, image):
        if chname != 'MOIRCS' or not self.gui_up:
            return

        imname = image.get('name', None)
        if imname is None:
            return

        # skip pre-baked "quicklook" (Q) frames
        if imname.startswith('MCSQ'):
            return

        header = image.get_header()

        frameid = header.get('FRAMEID', None)
        if frameid is None:
            return

        frameid = frameid.strip()
        fr = Frame(frameid)

        if 'DET-ID' not in header:
            return
        det_id = int(header['DET-ID'])
        # data_np = image.get_data()

        # if det_id == 1:
        #     # left chip?
        #     data_np = data_np[4:2044, 5:1821]
        #     hdu = fits.PrimaryHDU(data_np)
        #     img = AstroImage.AstroImage()
        #     img.load_hdu(hdu)
        #     img.set(name=f"{frameid}_QL", path=None, nothumb=True)
        #     channel = self.fv.get_channel_on_demand('MOIRCS_QL')
        #     channel.add_image(img)

        # elif det_id == 2:
        #     # right chip
        #     data_np = data_np[4:2044, 242:2010]
        #     hdu = fits.PrimaryHDU(data_np)
        #     img = AstroImage.AstroImage()
        #     img.load_hdu(hdu)
        #     img.set(name=f"{frameid}_QL", path=None, nothumb=True)
        #     channel = self.fv.get_channel_on_demand('MOIRCS_QL')
        #     channel.add_image(img)

        # TODO ...
        # can do some specialized processing on the image here

    def __str__(self):
        return 'ql_moircs'
