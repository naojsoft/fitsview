#
# FocusFit.py -- Focus fitting plugin for fits viewer
# 
# Yasu Sakakibara
# Eric Jeschke (eric@naoj.org)
#
import os.path

import gtk
import pango

import numpy

import gtk, gobject
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.figure as figure
from  matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as figurecanvas
from matplotlib.patches import Ellipse

import astro.fitsutils as fitsutils
import astro.curvefit as curvefit

from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
from ginga.gtkw import Plot
from ginga import GingaPlugin

class FocusFit(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(FocusFit, self).__init__(fv, fitsimage)

        self.lsf = curvefit.LeastSquareFits(self.logger)

    def build_gui(self, container):
        # Paned container is just to provide a way to size the graph
        # to a reasonable size
        box = gtk.VPaned()
        cw = container.get_widget()
        cw.pack_start(box, expand=True, fill=True)
        
        self.fontdesc = pango.FontDescription("Helvetica Bold 18")

        # Make the focus fitting plot
        self.fig = figure.Figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title('Focus Fitting')
        self.canvas = figurecanvas(self.fig)
        self.canvas.set_size_request(-1, 500)       
        #self.root.set_title("Focus Fitting")
        self.canvas.show_all()

        box.pack1(self.canvas, resize=True, shrink=True)

        # create a box to pack widgets into. 
        vbox1 = gtk.VBox()

        # label for seeing size 
        self.label_ss = gtk.Label("Seeing size: ")
        self.label_ss.set_alignment(0, 0)
        vbox1.pack_start(self.label_ss, padding=5, fill=True, expand=False)

        # label for data points
        self.label_dp = gtk.Label("Data points: ")
        self.label_dp.set_alignment(0, 0)
        vbox1.pack_start(self.label_dp, padding=5, fill=True, expand=False)
        #vbox1.show_all()

        fr = gtk.Frame(" QDAS Seeing ")
        fr.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
        fr.set_label_align(0.1, 0.5)
        fr.add(vbox1)
        fr.show_all()
        
        box.pack2(fr, resize=True, shrink=True)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(3)
        btns.set_child_size(15, -1)

        btn = gtk.Button("Close")
        btn.connect('clicked', lambda w: self.close())
        btns.add(btn)
        cw.pack_start(btns, padding=4, fill=True, expand=False)

    def close(self):
        chname = self.fv.get_channelName(self.fitsimage)
        self.fv.stop_local_plugin(chname, str(self))
        return True
        
    def clear(self):
        self.logger.debug('clearing canvas...')
        self.ax.cla()

    def _draw(self):
        self.ax.grid(True)
        self.fig.canvas.draw()

    def set_err_msg(self,msg,x,y):
        self.ax.text(x, y, msg, bbox=dict(facecolor='red',  alpha=0.1, ),
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=20, color='red')
  
    def display_error(self, msg):
        self.ax.annotate(msg, horizontalalignment='right', verticalalignment='bottom', fontsize=20)
        self._draw()


    def _drawGraph(self, title, result, data, minX, minY,
                   a, b, c):

        el = Ellipse((2, -1), 0.5, 0.5)

        self.clear()
        #self.set_title(title)
        self.ax.set_title('Focus Fitting')
        
        if result == 'empty':
            # all fwhm calculations failed; no data points found
            msg='No data available: FWHM all failed'
            self.set_err_msg(msg, 0.5, 0.5)
            self._draw()
            return False

        # <-- there are data points
        dpX  = [aDataPoint[0] for aDataPoint in data]
        dpY  = [aDataPoint[1] for aDataPoint in data]
        self.logger.debug("datapoints<%s> dpX<%s> dpY<%s>>" %(str(data), str(dpX),str(dpY) ))

        if result == 'singular':
            msg='bad data: result matrix is singular \ncurve fitting failed'
            self.ax.plot(dpX, dpY, 'ro')
            self.set_err_msg(msg, 0.5*(min(dpX)+max(dpX)), 0.5*(min(dpY)+max(dpY)))
            self._draw()
            return False

        if result == 'positive':
            msg="fitting curve calc failed\n a*x**2 + b*x + c \n a[%.3f] must be positive"  %(a)
            self.set_err_msg(msg, 0.465*(min(dpX)+max(dpX)), 0.5*(min(dpY)+max(dpY)))
            self.ax.plot(dpX, dpY, 'o')
            self._draw()
            return False

        # <-- maybe ok
        stX=dpX[0]
        enX=dpX[-1]
        numdataX = len(dpX)
        self.logger.debug(" stX[%s] minX[%s] enX[%s]  dpYmin<%s> minY<%s>  dpYmax<%s>" \
                           %(str(stX), str(minX), str(enX),  str(min(dpY)), str(minY), str(max(dpY)) ) )
                         
        # TODO: Hard coded constants!!
        t1 = numpy.arange(stX -0.035  , enX + 0.035, 0.005)
        
        t2 = self.lsf.plotQuadratic(t1, a, b, c)
        self.ax.plot(dpX, dpY, 'go', t1, t2)
                
        if dpX[0] <= minX <= dpX[-1]:
            self.ax.annotate('Z : %+.4f \nFWHM : %.4f(arcsec)' % (minX, minY) , xy=(minX, minY), xytext=(minX + 0.00001, max(dpY)-0.1),
                           bbox=dict(boxstyle="round",facecolor='green',ec="none",  alpha=0.1, ),
                           size=25, color='b',
                           arrowprops=dict(arrowstyle="wedge,tail_width=2.0",facecolor='green', ec="none", alpha=0.1, patchA=None, relpos=(0.5, -0.09)),
                           horizontalalignment='center')

               
            self.logger.debug("Z[%f] Fwhm[%f]" %(minX, minY) )
            #self.fig.canvas.set_window_title("Z[%f] Fwhm[%f]" %(minX, minY) )           

        else:
            self.logger.error("X<%s>  or Y<%s>  is out of range X<%s> Y<%s> " %(str(minX),str(minY), str(dpX), str(dpY)) )
            msg="z is out of range\nmin[%.3f] < z[%.3f] < max[%.3f]" %(dpX[0], minX, dpX[-1])
            self.set_err_msg(msg , 0.5*(min(dpX)+max(dpX)), 0.5*(min(dpY)+max(dpY)) )

        self._draw()
        return False

    def focus_fitting(self, file_list, x1, y1, x2, y2):
        try:
            # get beginning and ending frameids for title
            path, s_fits = os.path.split(file_list[0])
            path, e_fits = os.path.split(file_list[-1])
            s_fitsid, ext = os.path.splitext(s_fits)
            e_fitsid, ext = os.path.splitext(e_fits)
            title='%s ~ %s' % (s_fitsid, e_fitsid)
            self.logger.debug('fits %s' %(title)) 

        except OSError,e:
            self.logger.error('fail to set title %s' % str(e))
            title = ''
 
        z = None
        try:
            data_points = self.lsf.buildDataPoints(file_list, x1, y1, x2, y2)

            result = 'unknown'
            lsf_b = self.lsf.fitCurve(data_points)
            result = lsf_b.code

            z = lsf_b.minX; fwhm = lsf_b.minY
            self.logger.debug("result=%s z=%s fwhm=%s" % (result, z, fwhm))
            
            # draw graph at next available opportunity
            gobject.idle_add(self._drawGraph, title, result,
                             data_points, z, fwhm, lsf_b.a, lsf_b.b, lsf_b.c)

        except Exception, e:
            self.logger.error("focus fitting error <%s>" % (str(e)))

        return z           

    def seeing_size(self, avg, std, dp):
        if dp == 0:
            calc = False
            dp = 'Error'
        else:
            calc=True

        ss = "Seeing size: %5.2f+/-%5.2f(arcsec)" % (avg, std)
        self.label_ss.set_text(ss)
        
        dps = "Data points: %s" % str(dp)
        self.label_dp.set_text(dps)
        
        # setup font family, size, ....
        attr = pango.AttrList()
        if not calc:
            fg_color = pango.AttrForeground(65535, 0, 0, 0, 100)
            attr.insert(fg_color)
        else:
            fg_color = pango.AttrForeground(0, 0, 0, 0, 100)
            attr.insert(fg_color)
        self.label_ss.set_attributes(attr)
        self.label_dp.set_attributes(attr)  

        self.label_ss.modify_font(self.fontdesc)
        self.label_dp.modify_font(self.fontdesc)
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
        return 'focusfit'
    
# END
