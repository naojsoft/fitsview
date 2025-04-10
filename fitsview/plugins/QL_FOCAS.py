#
# QL_FOCAS.py -- QL_FOCAS plugin for Ginga reference viewer
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import os
import pathlib
import threading

from ginga import GingaPlugin
from ginga.misc import Bunch
from ginga.gw import Widgets, Viewers
from ginga.AstroImage import AstroImage

try:
    from naoj.focas import biassub
except ImportError:
    raise ImportError("Please install naojutils with the 'focas' bits")

from g2base.astro.frame import Frame

__all__ = ['QL_FOCAS']


import ObsLog

class QL_FOCAS(ObsLog.ObsLog):
    """
    QL_FOCAS
    ========
    FOCAS Gen2 quick look plugin.
    """

    def __init__(self, fv):
        super().__init__(fv)

        self.chnames = ['FOCAS_1', 'FOCAS_2']
        self.file_prefixes = ['FCSA']
        #self.sort_hdr = 'ExpID'

        # columns to be shown in the table
        column_info = [dict(col_title=self.sort_hdr, fits_kwd='FRAMEID'),
                       dict(col_title="Obs Mod", fits_kwd='OBS-MOD'),
                       dict(col_title="Datatype", fits_kwd='DATA-TYP'),
                       dict(col_title="Object", fits_kwd='OBJECT'),
                       dict(col_title="Date(UT)", fits_kwd='DATE-OBS'),
                       dict(col_title="Time(HST)", fits_kwd='HST-STR'),
                       #dict(col_title="PropId", fits_kwd='PROP-ID'),
                       dict(col_title="Exp Time", fits_kwd='EXPTIME'),
                       dict(col_title="Disperser", fits_kwd='DISPERSR'),
                       dict(col_title="Slit", fits_kwd='SLIT'),
                       dict(col_title="Bin Spatial", fits_kwd='BIN-FCT1'),
                       dict(col_title="Bin Spectrum", fits_kwd='BIN-FCT2'),
                       dict(col_title="Filter01", fits_kwd='FILTER01'),
                       dict(col_title="Filter02", fits_kwd='FILTER02'),
                       dict(col_title="Filter03", fits_kwd='FILTER03'),
                       dict(col_title="Air Mass", fits_kwd='AIRMASS'),
                       #dict(col_title="UT", fits_kwd='UT'),
                       dict(col_title="Pos Ang", fits_kwd='INST-PA'),
                       #dict(col_title="Ins Rot", fits_kwd='INSROT'),
                       #dict(col_title="Foc Val", fits_kwd='FOC-VAL'),
                       #dict(col_title="RA", fits_kwd='RA'),
                       #dict(col_title="DEC", fits_kwd='DEC'),
                       #dict(col_title="EQUINOX", fits_kwd='EQUINOX'),
                       dict(col_title="Memo", fits_kwd='G_MEMO'),
                       ]

        # Load preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QL_FOCAS')
        self.settings.set_defaults(sortable=True,
                                   color_alternate_rows=True,
                                   column_info=column_info)
        self.settings.load(onError='silent')

        self.col_info = self.settings.get('column_info', [])
        # this will set rpt_columns and col_widths
        self.process_columns(self.col_info)

        # construct path to where we are going to cache our quick look
        # result
        imdir = pathlib.Path(os.environ['GEN2COMMON']) / 'data_cache' / 'fitsview' / 'FOCAS'
        if not imdir.is_dir():
            imdir = None
        self.cache_dir = imdir

        # keeps track of exposures and their files
        self.exp_dct = dict()
        self.lock = threading.RLock()

    def build_gui(self, container):
        super().build_gui(container)

        proc_dir = os.path.join(os.environ['HOME'], 'Procedure')
        proc_dir_focas = os.path.join(proc_dir, 'FOCAS')
        if os.path.isdir(proc_dir_focas):
            proc_dir = proc_dir_focas
        self.w.obslog_dir.set_text(proc_dir)

        self.w.auto_save.set_state(True)

    def replace_kwds(self, header):
        d = super().replace_kwds(header)

        #d['FRAMEID'] = d['EXP-ID']
        return d

    def process_image(self, chname, header, image):
        if chname not in self.chnames:
            return

        imname = image.get('name', None)
        if imname is None:
            return

        # skip pre-baked "quicklook" (Q) frames
        if imname.startswith('FCSQ') or imname.endswith("_QL"):
            return

        path = image.get('path', None)
        if path is None:
            return

        header = image.get_header()
        exp_id, ch1_fits, ch2_fits = self.get_frames(path, header)
        if None in [exp_id, ch1_fits, ch2_fits]:
            # one of the frames has not arrived yet
            return

        if not self.gui_up:
            return
        self.reduce_ql(exp_id, ch1_fits, ch2_fits)

    def add_to_obslog(self, header, image):
        # if int(header.get('DET-ID', '')) != 1:
        #     return

        super().add_to_obslog(header, image)

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def start(self):
        #self.redo()
        pass

    ## def redo(self):
    ##     #self.bias_subtract_cb()
    ##     pass

    def stop(self):
        self.gui_up = False

    def reduce_ql(self, imname, ch1_fits, ch2_fits):
        # create a new image
        new_img = AstroImage(logger=self.logger)
        if self.cache_dir is None:
            impath = None
        else:
            impath = self.cache_dir / (imname + '.fits')
            # check if we have reduced this before--if so, just load
            # up our cached version
            if impath.exists():
                new_img.load_file(str(impath))
                new_img.set(name=imname)
                self.fv.gui_do(self.display_image, new_img)
                return

        try:
            hdulist = biassub.biassub(ch1_fits, ch2_fits)
            new_img.load_hdu(hdulist[0])

        except Exception as e:
            self.fv.show_error("Bias subtraction failed: %s" % (str(e)))
            return

        if impath is not None and not impath.exists():
            try:
                new_img.save_as_file(impath)
                new_img.set(path=str(impath), nothumb=False)
            except Exception as e:
                self.logger.warning(f"couldn't save {imname} as {impath}: {e}")
                new_img.set(path=None, nothumb=True)
        else:
            new_img.set(path=None, nothumb=True)

        new_img.set(name=imname)

        self.fv.gui_do(self.display_image, new_img)

    def get_frames(self, path, header):
        with self.lock:
            exp_id = header.get('EXP-ID', None)
            det_id = header.get('DET-ID', None)
            frame_id = header.get('FRAMEID', None)
            if None in [exp_id, det_id, frame_id]:
                return (exp_id, None, None)

            exp_id = exp_id.strip()
            det_id = str(int(det_id))
            frame_id = frame_id.strip()

            exp_bnch = self.exp_dct.get(exp_id, None)
            if exp_bnch is None:
                exp_bnch = Bunch.Bunch()
                self.exp_dct[exp_id] = exp_bnch

            det_bnch = exp_bnch.get(det_id, None)
            if det_bnch is None:
                det_bnch = Bunch.Bunch(path=path, frame_id=frame_id)
                exp_bnch[det_id] = det_bnch

            det1_bnch = exp_bnch.get('1', None)
            ch1_fits = None if det1_bnch is None else det1_bnch.path
            det2_bnch = exp_bnch.get('2', None)
            ch2_fits = None if det2_bnch is None else det2_bnch.path
            return (exp_id, ch1_fits, ch2_fits)

    def display_image(self, new_img):
        ch = self.fv.get_channel_on_demand("FOCAS_QL")
        ch.add_image(new_img)

    def __str__(self):
        return 'ql_focas'
