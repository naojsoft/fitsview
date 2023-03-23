#
# FOCAS.py -- FOCAS plugin for Ginga FITS viewer
#
# E. Jeschke
#
from ginga.misc import Widgets
from ginga.util import dp

from fitsview.util import focas
from fitsview.plugins import SPCAM


class FOCAS(SPCAM.SPCAM):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(FOCAS, self).__init__(fv, fitsimage)

        # Set preferences for destination channel
        prefs = self.fv.get_preferences()
        self.settings = prefs.createCategory('plugin_FOCAS')
        self.settings.setDefaults(annotate_images=False, fov_deg=0.15,
                                  match_bg=False, trim_px=0,
                                  merge=False, num_threads=4,
                                  drop_creates_new_mosaic=True,
                                  use_flats=False, flat_dir='',
                                  mosaic_new=True, make_thumbs=False)
        self.settings.load(onError='silent')

        self.dr = focas.FocasDR(logger=self.logger)

        self.mosaic_chname = 'FOCAS_Online'


    def get_exp_num(self, frame):
        exp_num = frame.number
        if (exp_num % 2) == 0:
            exp_num -= 1
        return exp_num

    def __str__(self):
        return 'focas'


#END
