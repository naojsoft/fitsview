#!/usr/bin/env python
#
# Takeshi Inagaki (tinagaki@naoj.org)
#
import os
import sys

import pygtk
pygtk.require('2.0')
import gtk
from gtk import gdk

import matplotlib
matplotlib.use('GTKAgg')  # or 'GTK'
# better rendering
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
# bad rendering #from matplotlib.backends.backend_gtk import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.figure import Figure

#from matplotlib import pylab

import numpy as np

import ssdlog


class CurveFitError(Exception):
    pass

class VertexError(CurveFitError):
    pass


class CurveFittingCanvas(FigureCanvas):

    def __init__(self, logger):

        self.fig = Figure(figsize=(5, 5), )

        super(CurveFittingCanvas, self).__init__(self.fig)
        self.axes = self.fig.add_subplot(111)
        self.set_size_request(-1, 400) 

        self.logger = logger


    def set_axes(self): 

        self.axes.grid(True)

        self.axes.set_xlabel('X')
        self.axes.set_ylabel('Y')

   
    def redraw(self):
        self.fig.canvas.draw()

    def clear_canvas(self):
        self.axes.cla()

    def error_message(self, msg, x_points, y_points):

        x = x_points.mean()
        y = y_points.mean()
        bbox_args = dict(boxstyle="round", fc="red", alpha=0.6)
        self.axes.annotate(msg, xy=(x, y), xytext=(x, y), 
                           size=20, color='g',
                           bbox=bbox_args,
                           ha='center',
                           va='bottom')

    def vertex(self, x, y):

        bbox_args = dict(boxstyle="round", fc="cyan", alpha=0.1)
        self.axes.annotate('Vertex(%.2f, %.2f)'%(x,y)  , xy=(x, y), xytext=(x, y), 
                           size=20, color='g',
                           bbox=bbox_args,
                           ha='center',
                           va='top')

    def plot_points(self, x_points, y_points):

        #self.axes.plot(x_points, y_points, 'b+', ms=10, mew=1.5)
        yerr = 0.005 * np.sqrt(x_points)
        self.axes.errorbar(x_points, y_points, yerr=yerr, fmt='s')

        #markerline, stemlines, baseline = pylab.stem(x_points, y_points, '-.')
        #pylab.setp(markerline, 'markerfacecolor', 'b')

        

        min_x = x_points.min()
        max_x = x_points.max()
        num = len(x_points)        

        margin = (max_x - min_x) / num
        
        self.axes.set_xlim([x_points[0]-margin, x_points[-1]+margin])

    def plot_graph(self, x_graph, y_graph):

        self.axes.plot(x_graph, y_graph, 'g-', linewidth=2)
        #self.canvas.draw()



class CurveFitting(CurveFittingCanvas):
    def __init__(self, qf, logger=None):
        super(CurveFitting, self).__init__(logger) 

        self.num_points = 10
        self.qf = qf
       
        self.logger = logger

    def parabola_downward(self):

        self.logger.debug('parabola opens downward...')

        try:
            max_x, max_y = self.qf.max_vertex()
            self.logger.debug('vertex(%f, %f)' %(max_x, max_y))

        except Exception as e:
            msg = 'error: coefficient A >= 0'
            raise VertexError(msg)

        else:
            return (max_x, max_y)        
            #self.error_message(msg, x, y)
            #self.logger.error('error: coefficient A is equal or greater than 0.')


    def parabola_upward(self):

        self.logger.debug('parabola opens upward...')

        try:
            min_x, min_y = self.qf.min_vertex()
            self.logger.debug('vertex(%f, %f)' %(min_x, min_y))

        except Exception as e:
            msg = 'error: coefficient A <= 0'
            raise VertexError(msg)
            #self.error_message(msg, x, y)
            #self.logger.error('error: coefficient A is equal or less than 0.')

        else:
            return (min_x, min_y)


    def plot(self, x_points, y_points, parabola, degree=2):

        self.logger.debug('start plotting...')
        #self.clear_canvas()

        self.set_axes()
        self.plot_points(x_points, y_points)

        opens = {'upward': self.parabola_upward, 
                 'downward': self.parabola_downward}

        try:
            self.qf.coefficient(x_points=x_points, y_points=y_points,
                                degree=degree)
            func = self.qf.quadratic()

        except Exception as e:
            msg = "error: failed to find quadratic equation"
            self.error_message(msg, x_points, y_points)
            errmsg = "error: failed to find quadratic equation: %s" % (str(e))
            self.logger.error(errmsg)
            raise CurveFitError(errmsg)   

        else:
           
            try:
                x, y = opens[parabola.lower()]()
                #max_x, max_y = self.qf.max_vertex()
                self.logger.debug('vertex(%f, %f)' %(x, y))

            except Exception as e:
                #msg = 'error: coefficient A >= 0'
                self.error_message(str(e), x_points, y_points)
                errmsg = 'error selecting parabola type: %s' % (str(e))
                self.logger.error(errmsg)
                raise CurveFitError(errmsg)   

            else:
                xs = len(x_points) * self.num_points
                x_new = np.linspace(x_points[0], x_points[-1], xs)
                y_new = func(x_new)
                #self.plot_points(x, y)
                self.plot_graph(x_new, y_new)
                self.vertex(x, y)
                return (x, y)

        #self.redraw()

def main(options,args):

    logname = 'curve_fitting'
    # Create top level logger.
    logger = ssdlog.make_logger(logname, options)

    # test data points
    xs = [5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600]

    if options.parabola == 'downward':
        # parabola opens downward
        ys = [11.598, 13.533, 13.564, 14.148, 13.443, 12.381, 10.997, 10.809, 8.695]
    
    if options.parabola == 'upward':
        # parabola opens upward
        ys = [11.598, 9.533, 7.564, 5.148, 3.443, 4.381, 6.997, 8.809, 10.695]
    
    # error case
    #ys = [11.598, 9.533, 23.564, 4.148, 0.443, 15.381, 10.997, 10.809, 20.695]
    #ys = [10, 10, 10, 10, 10, 10, 10, 10, 10]


    xpoints = np.asarray(xs)
    ypoints = np.asarray(ys)

    class MainWindow(gtk.Window):
 
        def __init__(self, logger):
            super(MainWindow, self).__init__()

            self.connect("destroy", lambda x: gtk.main_quit())
            self.set_default_size(600, 600)
            self.set_title("AO188 Curve Fitting!")

            from polynomial import QuadraticFunction
            self.qf = QuadraticFunction(self.logger)
            self.cf = CurveFitting(self.qf, logger=logger)

            self.set_gui()

        def set_gui(self): 

            vbox = gtk.VBox()
            self.add(vbox)
            vbox.pack_start(self.cf)

        def plot(self, x, y, parabola):
            self.cf.plot(x, y, parabola)

        def start(self):
            self.show_all()
            #pylab.show()
            gtk.main()

        def close(self):
            gtk.main_quit()


    # points = np.array([(1, -1), (1, 0.5), (1, 0.1), \
    #                    (2, 2), (2, 1.85), (2, 2.01), \
    #                    (3, 4.5), (3, 5.0), (3, 5.52), \
    #                    (4, 8), (4, 7.85), (4, 8.11), \
    #                    (5, 7), (5, 10.05), (5, 9.11), \
    #                    (6, 8.2), (6, 7.95), (6, 6.01), \
    #                    (7, 4.9), (7, 5.0), (7, 6.02), \
    #                    (8, 2), (8, 1.95), (8, 1.31), \
    #                    (9, 1.9), (9, 2.75), (9, 2.001),])
    # # get x and y vectors
    # x = points[:,0]
    # y = points[:,1]

    parabola = options.parabola

    try:
        main = MainWindow(logger=logger)
        main.plot(xpoints, ypoints, parabola)
        main.start()
    except KeyboardInterrupt:
        print 'keyboard interrupting...' 
        main.close()
        

if __name__ == "__main__":

    # Parse command line options with nifty new optparse module
    from optparse import OptionParser

    usage = "usage: %prog [options] [file] ..."
    optprs = OptionParser(usage=usage, version=('%prog'))

    optprs.add_option("-p", "--parabola", dest="parabola", default='downward',
                      help="Parabola opens upward|downward")
    
    optprs.add_option("--debug", dest="debug", default=False,
                      action="store_true",
                      help="Enter the pdb debugger on main()")

    optprs.add_option("--profile", dest="profile", action="store_true",
                      default=False,
                      help="Run the profiler on main()")

    ssdlog.addlogopts(optprs)
    (options, args) = optprs.parse_args(sys.argv[1:])

    #if len(args) > 0:
    #    optprs.error("incorrect number of arguments")
       
    # Are we debugging this?
    if options.debug:
        import pdb

        pdb.run('main(options, args)')

    # Are we profiling this?
    elif options.profile:
        import profile

        print "%s profile:" % sys.argv[0]
        profile.run('main(options, args)')

    else:
        main(options, args)






