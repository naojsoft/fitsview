# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
"""HSC Ginga QuickLook plugin.

**Plugin Type: Global**

``QL_HSC`` is a global plugin. Only one instance can be opened.

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
import threading


from naoj.hsc import hsc_dr
from g2base.astro.frame import Frame

import ObsLog

__all__ = ['QL_HSC']


class QL_HSC(ObsLog.ObsLog):

    def __init__(self, fv):
        super().__init__(fv)

        self.chname = 'HSC_Mosaic'
        self.frame_prefix = 'HSCA'
        self.data_dir = '/home/eric/testdata/HSC'

        self.dr = hsc_dr.HyperSuprimeCamDR(logger=self.logger)
        self.fov_deg = 2.0
        self.lock = threading.RLock()
        self.processed_frames = set([])
        self.sort_hdr = 'Exp_ID'
        self.sort_kwd = 'EXP-ID'

        # columns to be shown in the table
        column_info = [dict(col_title="Obs Mod", fits_kwd='OBS-MOD'),
                       dict(col_title=self.sort_hdr, fits_kwd=self.sort_kwd),
                       dict(col_title="Object", fits_kwd='OBJECT'),
                       dict(col_title="PropId", fits_kwd='PROP-ID'),
                       dict(col_title="HST", fits_kwd='HST-STR'),
                       #dict(col_title="OBCP Name", fits_kwd='T_UFNAME'),
                       #dict(col_title="UT", fits_kwd='UT'),
                       dict(col_title="Filter01", fits_kwd='FILTER01'),
                       dict(col_title="Exp Time", fits_kwd='EXPTIME'),
                       dict(col_title="Altitude", fits_kwd='ALTITUDE'),
                       dict(col_title="Azimuth", fits_kwd='AZIMUTH'),
                       dict(col_title="Air Mass", fits_kwd='AIRMASS'),
                       dict(col_title="Pos Ang", fits_kwd='INST-PA'),
                       dict(col_title="FocusZ", fits_kwd='FOC-VAL'),
                       #dict(col_title="Ins Rot", fits_kwd='INSROT'),
                       dict(col_title="RA", fits_kwd='RA'),
                       dict(col_title="DEC", fits_kwd='DEC'),
                       dict(col_title="EQUINOX", fits_kwd='EQUINOX'),
                       dict(col_title="Memo", fits_kwd='G_MEMO'),
                       ]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QL_HSC')
        self.settings.set(sortable=True,
                          color_alternate_rows=True,
                          column_info=column_info,
                          cache_normalized_images=True)
        self.settings.load(onError='silent')

        self.col_info = self.settings.get('column_info', column_info)
        # this will set rpt_columns and col_widths
        self.process_columns(self.col_info)

    def replace_kwds(self, dct):
        exp_id = dct.get(self.sort_kwd, None)
        if exp_id is None:
            return dct

        if 'G_MEMO' not in dct:
            dct['G_MEMO'] = ''
        return dct

    def start(self):
        super().start()

        channel = self.fv.get_channel_on_demand('HSC_Mosaic')
        #channel.fitsimage.add_callback('drag-drop', self.drop_cb)

    def stop(self):
        super().stop()

    def incoming_data_cb(self, fv, chname, image, info):
        imname = image.get('name', None)
        if imname is None:
            return

        # only accepted list of frames
        if not imname.startswith(self.frame_prefix):
            return

        fr = Frame(image.get('path', imname))
        exp_num = self.dr.get_exp_num(fr.frameid)

        with self.lock:
            # have we processed this image before
            if exp_num in self.processed_frames:
                return
            self.processed_frames.add(exp_num)

        self.logger.info(f"incoming exposure {exp_num}")
        header = image.get_header()
        dct = self.replace_kwds(header)

        # add image to obslog
        self.fv.gui_do(self.add_to_obslog, dct, image)

    # def process_image(self, chname, header, image):
    #     if chname != 'HSC':
    #         return

    #     imname = image.get('name', None)
    #     if imname is None:
    #         return

    #     header = image.get_header()

    #     frameid = header.get('FRAMEID', None)
    #     if frameid is None:
    #         return

    #     frameid = frameid.strip()

    #     # normalized image prefix
    #     newname = 'HSCE' + frameid[4:-2] + '00'


    def preprocess(self, image):
        # TODO: have this remove overscan regions, and anything else
        # desirable
        self.dr.remove_overscan(image, sub_bias=True)
        return image

    def build_collage(self, exp_num, data_dir=None):
        self.logger.info(f"building collage for exp id '{exp_num}'")
        if data_dir is None:
            data_dir = self.data_dir
        channel = self.fv.get_current_channel()
        if channel.name != self.chname:
            channel = self.fv.get_channel_on_demand(self.chname)
            self.fv.change_channel(self.chname)

        # start plugin on HSC_Mosaic channel
        plugin_name = 'Collage'
        opmon = channel.opmon
        if not opmon.is_active(plugin_name):
            self.fv.start_local_plugin(self.chname, plugin_name, None)

        paths = self.dr.exp_num_to_file_list(data_dir, exp_num)

        plugin = opmon.get_plugin(plugin_name)
        plugin.collage(paths, preprocess=self.preprocess, new_collage=True)

    def dblclick_cb(self, widget, d):
        """Collage the image that was double-clicked in the obslog"""
        key = list(d.keys())[0]
        exp_num = int(key[4:] + '00')
        self.logger.info(f"User clicked on exposure {exp_num}")

        self.build_collage(exp_num)
        return True

    def __str__(self):
        return 'ql_hsc'
