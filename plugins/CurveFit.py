#
# CurveFit.py -- Curve fitting plugin for fits viewer
# 
# Takeshi Inagaki (tinagaki@naoj.org)
# Eric Jeschke (eric@naoj.org)
#
import os.path

#import gtk
#import pango

import numpy

#import gtk

import Gen2.fitsview.util.curve_fit as curvefit
#from Gen2.fitsview.util.polynomial import QuadraticFunction

from ginga.misc import CanvasTypes, Widgets
from ginga.gw import Plots

#from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
#from ginga.gtkw import Plot
from ginga import GingaPlugin

class CurveFit(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(CurveFit, self).__init__(fv, fitsimage)

        #self.qf = QuadraticFunction(self.logger)

    def build_gui(self, container, future=None):
        # Paned container is just to provide a way to size the graph
        # to a reasonable size
        self.logger.debug('building curve fitting...')

        self.plot1 = Plots.Plot(logger=self.logger,
                                width=2, height=6, dpi=72)

        class1 = curvefit.make_CurveFittingCanvas(self.plot1.fig.canvas.__class__)
        class2 = curvefit.make_CurveFitting(class1)
        self.cf = class2(self.plot1.fig, logger=self.logger)
        # GTK only
        self.cf.set_size_request(-1, 700)
        
        cf = Widgets.wrap(self.cf)
  
        #self.cf = cf.get_widget()
        #self.t.clear_canvas()

        vtop = Widgets.VBox()
        vtop.set_border_width(2)
        #box = gtk.VPaned()
        #cw = container.get_widget()
        #cw.pack_start(box, expand=True, fill=True)
        #cw.pack_start(self.cf, expand=True, fill=True)
        #vbox, sw, orientation = Widgets.get_oriented_box(container)
        box = Widgets.Splitter(orientation='vertical')
        
        #self.fontdesc = pango.FontDescription("Helvetica Bold 18")

        #box.pack1(self.cf, resize=True, shrink=True)
        #cfbox = Widgets.VBox()
        #cfbox.add_widget(self.cf)
        #box.add_widget(cfbox)
        box.add_widget(cf)

        #btns = gtk.HButtonBox()
        btns = Widgets.HBox()
        #btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(3)
        #btns.set_child_size(15, -1)

        # create an empty box to adjust height of cure fitting. 
        #vbox1 = gtk.VBox()
        vbox1 = Widgets.VBox()
        vbox1.set_border_width(100)

        # empty label
        #self.label = gtk.Label("")
        self.label = Widgets.Label("")
        #self.label.set_alignment(0, 0)
        #vbox1.pack_start(self.label, padding=50, fill=True, expand=False)
        vbox1.add_widget(self.label)

        #fr = gtk.Frame()
        fr = Widgets.Frame("")
        #fr.set_shadow_type(gtk.SHADOW_NONE)
        #fr.set_label_align(0.1, 0.5)
        #fr.add(vbox1)
        fr.set_widget(vbox1)
        #fr.show_all()

        #box.pack2(fr, resize=True, shrink=True)
        box.add_widget(fr)
        #vbox.add_widget(box, stretch=1)

        #btn = gtk.Button("Close")
        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        #btns.add(btn)
        btns.add_widget(btn, stretch=0)

        #cw.pack_start(btns, padding=4, fill=True, expand=False)
        #vtop.add_widget(sw, stretch=1)
        vtop.add_widget(box, stretch=1)
        vtop.add_widget(btns, stretch=0)
        container.add_widget(vtop, stretch=1)

        #container.add_widget(box)
        #container.add_widget(btns)

    def close(self):
        chname = self.fv.get_channelName(self.fitsimage)
        self.fv.stop_local_plugin(chname, str(self))
        return True
 
    def curve_fitting(self, p, x_points, y_points, parabola):

        x = numpy.asarray(x_points)
        y = numpy.asarray(y_points)

        try:
            self.cf.clear_canvas()
            self.cf.set_axes()
            mx, my = self.cf.plot(x, y, parabola)
            self.cf.redraw()

            p.setvals(result='ok', mx=float(mx), my=float(my),
                      a=float(self.cf.qf.a), b=float(self.cf.qf.b),
                      c=float(self.cf.qf.c))

        except Exception as e:
            errmsg = "error in curve fitting: %s" % (str(e))
            self.logger.error(errmsg)
            p.setvals(result='error', errmsg=str(e))

        return 0

    def start(self):
        self.resume()

    def pause(self):
        pass
        
    def resume(self):
        pass
        
    def stop(self):
        self.fv.showStatus("")
        
    def redo(self):
        pass
    
    def __str__(self):
        return 'curvefit'
    
# END
