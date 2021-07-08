#
# QL_FOCAS.py -- QL_FOCAS plugin for Ginga reference viewer
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
import os

from ginga import GingaPlugin
from ginga.gw import Widgets, Viewers
from ginga.AstroImage import AstroImage

try:
    from naoj.focas import biassub, reconstruct_image as recon
except ImportError:
    raise ImportError("Please install naojutils with the 'focas' bits")

from astro.frame import Frame

__all__ = ['QL_FOCAS']


class QL_FOCAS(GingaPlugin.GlobalPlugin):
    """
    QL_FOCAS
    ======
    FOCAS Gen2 quick look plugin.

    Usage
    -----
    """

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(QL_FOCAS, self).__init__(fv)

        self._wd = 300
        self._ht = 300
        self.q_image = None
        self.fitsimage = None

        # Load preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_QL_FOCAS')
        self.settings.set_defaults(region_file=None)
        self.settings.load(onError='silent')

        self.region_file = self.settings.get('region_file', None)

    def build_gui(self, container):
        # Users sometimes don't open the plugin on the correct channel.
        # Force the correct channel to be used.
        chinfo = self.fv.get_channel_on_demand('FOCAS')
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
        zi.enable_autozoom('override')
        zi.enable_autocuts('override')
        zi.set_autocenter('on')
        #zi.set_zoom_algorithm('step')
        zi.set_zoom_algorithm('rate')
        zi.set_zoomrate(1.62)
        zi.show_mode_indicator(True)
        zi.show_color_bar(True)
        settings = zi.get_settings()
        zi.set_bg(0.4, 0.4, 0.4)
        zi.set_color_map('gray')
        zi.set_color_map('gray')
        zi.set_intensity_map('ramp')
        # for debugging
        zi.set_name('focas_qimage')
        zi.add_callback('cursor-changed', self.motion_cb)
        self.q_image = zi

        # add a canvas for drawing the slices
        canvas = zi.get_canvas()
        DrawingCanvas = canvas.get_draw_class('drawingcanvas')
        self.canvas = DrawingCanvas()
        canvas.add(self.canvas)

        bd = zi.get_bindings()
        bd.enable_all(True)

        iw = Viewers.GingaViewerWidget(zi)
        iw.resize(self._wd, self._ht)
        vbox1.add_widget(iw, stretch=1)

        fr = Widgets.Frame("Reduced")
        fr.set_widget(vbox1)
        top.add_widget(fr, stretch=1)

        fr = Widgets.Frame("FOCAS")

        captions = (('Exposure:', 'label', 'Exposure', 'entry'),
                    ('Bias Subtract', 'button',
                     'Mark slices', 'button',
                     'IFU Quick Look', 'button',
                     'Insert Image', 'button'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w = b

        b.bias_subtract.add_callback('activated', self.bias_subtract_cb)
        b.bias_subtract.set_tooltip("Bias subtract and combine this exposure")
        b.mark_slices.add_callback('activated', self.mark_slices_cb)
        b.mark_slices.set_tooltip("Mark the slices on the image")
        b.ifu_quick_look.add_callback('activated', self.ifu_quick_look_cb)
        b.ifu_quick_look.set_tooltip("IFU quick look this exposure")
        b.insert_image.add_callback('activated', self.insert_image_cb)
        b.insert_image.set_tooltip("Add image to FOCAS_QL channel")

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

    def get_frames(self):
        exp = self.w.exposure.get_text().strip()
        if len(exp) == 0:
            image = self.fitsimage.get_image()

            if image is None:
                # Nothing to do
                self.fv.show_error("No image to quick reduce!")
                return

            path = image.get('path', None)
            if path is None:
                self.fv.show_error("No path in image!")
                return
        else:
            fr = Frame('FCSA%08d.fits' % int(exp))
            fr.directory = os.environ['DATAHOME']
            path = fr.path

        fr = Frame(path)
        exp_num = fr.number - (fr.number % 2)

        fr.number = exp_num
        ch1_fits = fr.path

        fr.number = exp_num + 1
        ch2_fits = fr.path
        return (ch1_fits, ch2_fits)

    def bias_subtract_cb(self, w):
        self.canvas.delete_all_objects()

        ch1_fits, ch2_fits = self.get_frames()
        fr = Frame(ch1_fits)
        name = str(fr) + '_QL'

        self.q_image.onscreen_message("Working ...")

        try:
            hdulist = biassub.bias_subtraction(ch1_fits, ch2_fits)

            # create a new image
            new_img = AstroImage(logger=self.logger)
            new_img.load_hdu(hdulist[0])

            # no thumbnails presently
            new_img.set(nothumb=True, path=None, name=name)

            self.q_image.set_image(new_img)

        except Exception as e:
            self.fv.show_error("Bias subtraction failed: %s" % (str(e)))

        finally:
            self.q_image.onscreen_message(None)

    def ifu_quick_look_cb(self, w):
        self.canvas.delete_all_objects()

        ch1_fits, ch2_fits = self.get_frames()
        fr = Frame(ch1_fits)
        name = str(fr) + '_IFU'

        self.q_image.onscreen_message("Working ...")

        try:
            reg_data = recon.read_region_file(self.region_file)

            hdulist = recon.reconstruct_image(ch1_fits, ch2_fits,
                                              reg_data)

            # create a new image
            new_img = AstroImage(logger=self.logger)
            new_img.load_hdu(hdulist[0])

            # no thumbnails presently
            new_img.set(nothumb=True, path=None, name=name)

            self.q_image.set_image(new_img)

        except Exception as e:
            self.fv.show_error("IFU quick look failed: %s" % (str(e)))

        finally:
            self.q_image.onscreen_message(None)


    def insert_image_cb(self, w):
        channel = self.fv.get_channel_on_demand('FOCAS_QL')
        img = self.q_image.get_image()
        channel.add_image(img)

    def mark_slices_cb(self, w):
        self.canvas.delete_all_objects()

        Box = self.canvas.get_draw_class('box')
        Text = self.canvas.get_draw_class('text')

        reg_data = recon.read_region_file(self.region_file)

        for i in range(len(reg_data)):
            x, y, wd, ht = reg_data[i]
            x, y, xr, yr = x - 1, y - 1, wd / 2.0, ht / 2.0
            self.canvas.add(Box(x, y, xr, yr, color='green'),
                            redraw=False)
            num = len(reg_data) - i
            self.canvas.add(Text(x, y+yr+6, str(num), color='green'),
                            redraw=False)
        self.canvas.update_canvas(whence=3)

    ## def redo(self):
    ##     #self.bias_subtract_cb()
    ##     pass

    def motion_cb(self, viewer, button, data_x, data_y):
        self.fv.showxy(viewer, data_x, data_y)
        return True

    def __str__(self):
        return 'ql_focas'


#END
