#
# QL_PFS.py -- PFS quick look plugin for fits viewer
#
# C. Loomis
# E. Jeschke
#
"""Do a quick look for PFS images.

**Plugin Type: Global**

``QL_PFS`` is a global plugin. Only one instance can be opened.

**Usage**

"""
import os
import pathlib

import numpy as np
from astropy.io import fits

#from ginga.gw import Widgets
from ginga import GingaPlugin
from ginga.AstroImage import AstroImage


class QL_PFS(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        super().__init__(fv)

        # prefs = self.fv.get_preferences()
        # self.settings = prefs.create_category('plugin_QL_PFS')
        # self.settings.add_defaults()
        # self.settings.load(onError='silent')

        # construct path to where we are going to cache our quick look
        # result
        imdir = pathlib.Path(os.environ['GEN2COMMON']) / 'data_cache' / 'PFS'
        if not imdir.is_dir():
            imdir = None
        self.cache_dir = imdir

        # No UI, at present
        #self.gui_up = False

    def start(self):
        pass

    def stop(self):
        #self.gui_up = False
        pass

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def redo(self, channel, image):
        # if not self.gui_up:
        #     return

        if image is None:
            return
        path = image.get('path', None)
        if path is None:
            self.logger.info("No path set for image--skipping")
            return

        name = image.get('name', None)
        if name is None:
            self.logger.info("No name set for image--skipping")
            return

        if not name.startswith('PFS'):
            self.logger.debug("Not a PFS file--skipping")
            return

        if name.startswith('PFSB') and not '_QL' in name:
            #a_img = self._reduce_ql(path)
            self.fv.nongui_do(self._reduce_ql, channel, path)
        else:
            self.logger.debug("Not a raw PFS 'B' file--nothing to do")

    def display_image(self, channel, a_img):
        channel.add_image(a_img)
        # TODO: create the channel in the same workspace
        # chname = channel.name + '_QL'
        # self.fv.gui_do(self.fv.add_image, a_img.get('name'), a_img,
        #                chname=chname)

    # spun off into a different function so we can run it in a different
    # thread or process if needed
    #
    # @cloomis's pseudo-code:
    # nreads = PHDU['W_H4NRED']
    # readN = hdulist[f'IMAGE_{nreads}'].astype('f4') - hdulist[f'REF_{nreads}']
    # read1 =  hdulist[f'IMAGE_1'].astype('f4') - hdulist[f'REF_1']
    # rampCdsImage = readN - read1
    #
    def _reduce_ql(self, channel, path):
        p = pathlib.Path(path)

        a_img = AstroImage(logger=self.logger)
        imname = f"{p.stem}_QL"
        if self.cache_dir is None:
            impath = None
        else:
            impath = self.cache_dir / (imname + '.fits')
            # check if we have reduced this before--if so, just load
            # up our cached version
            if impath.exists():
                a_img.load_file(str(impath))
                a_img.set(name=imname)
                self.fv.gui_do(self.display_image, channel, a_img)
                return

        a_img.set(name=imname, path=str(impath))

        # use astropy.io.fits because we can explicitly set memmap and
        # lazy loading of HDUs
        with fits.open(path, 'readonly', memmap=False,
                       lazy_load_hdus=True) as pfsb_f:
            nreads = pfsb_f[0].header.get('W_H4NRED', None)
            if nreads is None:
                # just show last HDU
                a_img.load_hdu(pfsb_f[-1])
            else:
                read_n = (pfsb_f[f'IMAGE_{nreads}'].data.astype('f4') -
                          pfsb_f[f'REF_{nreads}'].data)
                read_1 = (pfsb_f[f'IMAGE_1'].data.astype('f4') -
                          pfsb_f[f'REF_1'].data)
                ramp_img = read_n - read_1
                a_img.set_data(ramp_img)
                # what header should we use for this?
                a_img.update_keywords(pfsb_f[f'IMAGE_{nreads}'].header)

            # force close so no lingering
            try:
                pfsb_f.close()
            finally:
                pass

            if impath is not None and not impath.exists():
                try:
                    a_img.save_as_file(impath)
                except Exception as e:
                    self.logger.warning(f"couldn't save {imname} as {impath}: {e}")
                    a_img.set(path=None)

        self.fv.gui_do(self.display_image, channel, a_img)

    def __str__(self):
        return 'ql_pfs'
