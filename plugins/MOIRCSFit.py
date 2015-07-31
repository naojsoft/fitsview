#
# MOIRCSFit.py -- MOIRCS focus fitting plugin for fits viewer
# 
# Hannah Twigg-Smith
# Writen using FocusFit.py as a guideline
#

import os

import gtk
import pango

import numpy as np

import gtk, gobject
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
from  matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as figurecanvas
import pylab

from astropy.table import Table

import pickle
import sewpy
import pyfits

from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
from ginga.gtkw import Plot
from ginga import GingaPlugin

class MOIRCSFit(GingaPlugin.LocalPlugin):
    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(MOIRCSFit, self).__init__(fv, fitsimage)

    def build_gui(self, container):
        # Paned container is just to provide a way to size the graph
        # to a reasonable size
        box = gtk.VPaned()
        cw = container.get_widget()
        cw.pack_start(box, expand=True, fill=True)

        self.notebook = gtk.Notebook()
        self.notebook.set_tab_pos(gtk.POS_TOP)
        self.show_tabs = True
        self.show_border = True
        self.notebook.set_scrollable(True)
        box.pack1(self.notebook, resize=True, shrink=True)
        
        self.fontdesc = pango.FontDescription("Helvetica Bold 18")

        # create a box to pack widgets into. 
        vbox1 = gtk.VBox()

        # label for seeing size 
        self.label_ss = gtk.Label(" Best: ")
        self.label_ss.set_alignment(0, 0)
        vbox1.pack_start(self.label_ss, padding=5, fill=True, expand=False)

        fr = gtk.Frame(" Best Z Value ")
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
        self.ax[0].cla()
        self.ax[1].cla()

    def _draw(self):
        self.ax[0].grid(True)
        self.ax[1].grid(True)
        self.fig.canvas.draw()

    def _drawGraph(self, x_pixels, y_fwhm, fits_names, FOCZ):
        self.fig, self.ax = plt.subplots(2, 1)
        self.canvas = figurecanvas(self.fig)
        self.canvas.set_size_request(-1, 500)

        for i in range(0, 2):
            fit = np.polyfit(x_pixels[i*2+1],y_fwhm[i*2+1],1)
            fit_fn = np.poly1d(fit)

            self.ax[i].plot(x_pixels[i*2], y_fwhm[i*2], 'ro', label='Raw data')
            self.ax[i].plot(x_pixels[i*2+1], y_fwhm[i*2+1], 'go', label='Selected')
            self.ax[i].plot(x_pixels[i*2+1], fit_fn(x_pixels[i*2+1]), 'k', label='Best fit')
            self.ax[i].set_title('%s: FOC-Z=%s, Slope=%.3g' %(fits_names[i], FOCZ[i], fit[0]))
            self.ax[i].axis([0, 2048, 0, 2.5])
            self.ax[i].set_xlabel('x [pixel]')
            self.ax[i].set_ylabel('FWHM [arcsec]')
            self.ax[i].legend(fontsize='small')
            self.ax[i].set_xticks(np.arange(0, 2049, 512))

        label = gtk.Label('Z='+str(FOCZ[0]))
        self.notebook.append_page(self.canvas,label)
        self._draw()
        self.canvas.show_all()
        self.logger.info('%s: FOC-Z=%s, Slope=%.3g' %(fits_names[i], FOCZ[i], fit[0]))
        self.notebook.next_page()
        return False

    def errorplot(self, det_id, astro_data):

        x = []
        y = []
        e = []

        for i in astro_data:
            x.append(i[det_id-1][2])
            y.append(np.mean(i[det_id-1][4]))
            e.append(np.mean(i[det_id-1][4])-np.sqrt(np.mean(np.square(i[det_id-1][4])))) # RMS standard deviation

        z = np.polyfit(x,y,2)   # Fit a parabola to the xy points
        fit_fn = np.poly1d(z)
        best = z[1]/(-2.0*z[0]) # Find the y value of the vertex of the 
                                # parabola, (-b/2a). This is the best focus
                                # position for this detector.

        self.ax[det_id-1].errorbar(x, y, e, linestyle='None', marker='^')
        xpoints = np.linspace(x[0],x[-1],1000)
        self.ax[det_id-1].plot(xpoints, fit_fn(xpoints))

        self.ax[det_id-1].set_title('MCS Focus Fitting [CH%s]: Best-Z=%.3g' %(det_id, best))
        self.ax[det_id-1].set_xlabel('FOC-Z')
        self.ax[det_id-1].set_ylabel('FWHM [arcsec]')
        self.logger.info('MCS Focus Fitting [CH%s]: Best-Z=%.3g' %(det_id, best))
        return best

    def _drawBest(self):
        
        try:
            astro = pickle.load(open('/tmp/astro_data.pck','r'))
            self.logger.debug('Unpickled!')
        
        except Exception, e:
            self.logger.error('Unpickling failed')
            astro = []

        astro.sort()

        self.fig, self.ax = plt.subplots(2, 1)
        self.canvas = figurecanvas(self.fig)
        self.canvas.set_size_request(-1, 500)

        self.clear()

        bests = [self.errorplot(1, astro), self.errorplot(2, astro)]

        self.label_ss.set_text(" Best:  %.3g" %(np.average(bests)))
        self._draw()
        self.canvas.show_all()
        label = gtk.Label("Best")
        self.notebook.append_page(self.canvas,label)
        self.notebook.next_page()
        os.system('rm /tmp/astro_data.pck')

        return None

    def valuefilter(self, astrotable):
        indices = []
        for row in astrotable:
            row['FWHM_IMAGE'] *= 0.117 # Multiply the FWHM row by 0.117
        for row in astrotable:
            if (row['FLAGS'] >= 2 or
                row['ELLIPTICITY'] >= 0.3 or
                row['FWHM_IMAGE'] == 0):
                indices.append(row.index)
        astrotable.remove_rows(indices) # Throw out unneeded values
        return astrotable

    def adjust(self, astrotable):
        for row in astrotable:
            row['FWHM_IMAGE'] *= 0.117 # Multiply the FWHM row by 0.117
        return astrotable

    def focus_fitting(self, param, propid, path1, path2, fits1, fits2):
        fits_names = [fits1, fits2]
        path_names = [path1, path2]
        DETID = []
        FOCZ = []
        
        try:
            for i in path_names:
                hdulist = pyfits.open(i)
                DETID.append(hdulist[0].header['DET-ID']) # Get Detector ID #s
                FOCZ.append(hdulist[0].header['FOC-VAL']) # Get Focus Values
                hdulist.close()
            self.logger.debug('Got DETID and FOCZ from fits headers')

        except Exception, e:
             self.logger.error('Fail to pull DETID and FOCZ from fits headers')

        detector_pair = [] #Initialize list that will be pickled
        for i in range(0,len(fits_names)):
            detector_pair.append([fits_names[i],DETID[i],FOCZ[i]])

        sew = sewpy.SEW(
		params=["X_IMAGE", "Y_IMAGE", "FWHM_IMAGE", "FLAGS",
                        "ELLIPTICITY", "NUMBER"],
		configfilepath='/home/hannahts/Code/pymcsfcs/focus-CONFIG.sex',
		sexpath='/home/hannahts/Code/sextractor-2.8.6/src/sex',
                loglevel='ERROR'
         	)
        
        catalog = []
        
        try:
            for i in path_names:
                table = sew(i)['table']
                table2 = sew(i)['table']
                catalog.append(self.adjust(table))
                catalog.append(self.valuefilter(table2))
            self.logger.debug('Source Extractor')

        except Exception, e:
            self.logger.error('Source Extractor failed:')

        x_pixels=[]
        y_fwhm = []

        for table in catalog:
            int_x = []
            int_y = []
            for row in table:
                int_x.append(row['X_IMAGE'])
                int_y.append(row['FWHM_IMAGE'])
            x_pixels.append(int_x)
            y_fwhm.append(int_y)

        for i in range(0, 2):
            detector_pair[i].append(x_pixels[i*2+1])
            detector_pair[i].append(y_fwhm[i*2+1])

        try:
            if os.path.isfile('/tmp/astro_data.pck'):
                astro_data = pickle.load(open('/tmp/astro_data.pck','rb'))
                astro_data.append(detector_pair)
                f = open('/tmp/astro_data.pck','wb')
                pickle.dump(astro_data, f)
                f.close()
            else:
                astro_data = [detector_pair]
                f = open('/tmp/astro_data.pck','wb')
                pickle.dump(astro_data, f)
                f.close()
            self.logger.debug('pickling complete')
        
        except Exception, e:
            self.logger.error('pickling error')

        z = None
        try:
            # draw graph
            self._drawGraph(x_pixels, y_fwhm, fits_names, FOCZ)

        except Exception, e:
            self.logger.error("focus fitting error")

        return z

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
        return 'moircsfit'
    
# END
