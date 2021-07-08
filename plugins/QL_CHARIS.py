#
# QL_CHARIS.py -- QuickLook CHARIS plugin for Ginga reference viewer
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
from ginga import GingaPlugin
from ginga.gw import Widgets, Viewers
from ginga.AstroImage import AstroImage

import numpy as np
from astropy.io import fits

class QL_CHARIS(GingaPlugin.GlobalPlugin):
    """
    QL_CHARIS
    ======
    CHARIS Gen2 quick look plugin.

    Usage
    -----
    """

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(QL_CHARIS, self).__init__(fv)

        self.sb_hdu1 = True
        self.hdu1_path = None
        self.hdu1_data = None

        self._wd = 300
        self._ht = 300
        self.q_image = None
        self.fitsimage = None

    def build_gui(self, container):
        # Users sometimes don't open the plugin on the correct channel.
        # Force the correct channel to be used.
        chinfo = self.fv.get_channel_on_demand('CHARIS')
        self.fitsimage = chinfo.fitsimage

        top = Widgets.VBox()
        top.set_border_width(4)
        top.set_spacing(2)

        vbox1 = Widgets.VBox()

        # Uncomment to debug; passing parent logger generates too
        # much noise in the main logger
        #zi = Viewers.CanvasView(logger=self.logger)
        zi = Viewers.CanvasView(logger=None)
        zi.set_desired_size(self._wd, self._ht)
        zi.enable_autozoom('once')
        zi.enable_autocuts('once')
        #zi.enable_autocenter('once')
        zi.set_autocenter('once')
        zi.set_zoom_algorithm('step')
        zi.show_mode_indicator(True)
        zi.show_color_bar(True)
        settings = zi.get_settings()
        zi.set_bg(0.4, 0.4, 0.4)
        zi.set_color_map('gray')
        zi.set_color_map('gray')
        zi.set_intensity_map('ramp')
        # for debugging
        zi.set_name('charis_qimage')
        zi.add_callback('cursor-changed', self.motion_cb)
        self.q_image = zi

        bd = zi.get_bindings()
        bd.enable_all(True)

        iw = Viewers.GingaViewerWidget(zi)
        iw.resize(self._wd, self._ht)
        vbox1.add_widget(iw, stretch=1)

        fr = Widgets.Frame("Reduced")
        fr.set_widget(vbox1)
        top.add_widget(fr, stretch=1)

        fr = Widgets.Frame("Charis")

        captions = (('Subtract HDU 0', 'checkbutton'),
                    ('Quick Reduce', 'button'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w = b

        chk_btn = b.subtract_hdu_0
        chk_btn.set_state(self.sb_hdu1)
        chk_btn.add_callback('activated', self.toggle_sb_hdu1)

        b.quick_reduce.add_callback('activated', lambda w: self.quick_reduce())
        b.quick_reduce.set_tooltip("Update from the current image in the channel")

        fr.set_widget(w)
        top.add_widget(fr, stretch=0)

        #spacer = Widgets.Label('')
        #top.add_widget(spacer, stretch=1)

        btns = Widgets.HBox()
        btns.set_spacing(3)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        ## btn = Widgets.Button("Help")
        ## btn.add_callback('activated', lambda w: self.help())
        ## btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)

    def toggle_sb_hdu1(self, w, tf):
        self.sb_hdu1 = tf

        self.redo()

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def start(self):
        #self.redo()
        pass

    def pause(self):
        pass

    def resume(self):
        pass

    def stop(self):
        self.pause()
        self.hdu1_path = None
        self.hdu1_data = None

    def quick_reduce(self):
        image = self.fitsimage.get_image()
        if image is None:
            # Nothing to do
            return

        path = image.get('path', None)
        if path is None:
            return

        self.q_image.onscreen_message("Working ...")

        try:
            #img_data = image.get_data().astype(np.float32)
            with fits.open(path, 'readonly') as in_f:
                n = len(in_f) - 1
                img_data = in_f[n].data.astype(np.float32)

                if self.sb_hdu1:
                    # get HDU 1 from the file and subtract it

                    if path != self.hdu1_path:
                        # if not cached then we have to re-fetch it
                        self.hdu1_path = path
                        self.hdu1_data = in_f[1].data

                    sbr_data = img_data - self.hdu1_data

                else:
                    sbr_data = img_data

                # create a new image
                metadata = dict(header=image.get_header())
                new_img = AstroImage(data_np=sbr_data, metadata=metadata,
                                     logger=self.logger)
                # no thumbnails presently
                new_img.set(nothumb=True, path=None, name=image.get('name'))

                self.q_image.set_image(new_img)

        finally:
            self.q_image.onscreen_message(None)


    ## def redo(self):
    ##     #self.quick_reduce()
    ##     pass

    def motion_cb(self, viewer, button, data_x, data_y):
        self.fv.showxy(viewer, data_x, data_y)
        return True

    def __str__(self):
        return 'ql_charis'


#END
