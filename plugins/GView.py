#
# GView.py -- plugin for Ginga implementing some of the commands from
#               the old ZView viewer
#
# Eric Jeschke (eric@naoj.org)
#
# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
#
"""

Examples:
    rd 1 /home/eric/testdata/SPCAM/SUPA01118760.fits
    v 1
    cm jt
    cd /home/eric/testdata/SPCAM
    rd 1 SUPA01118760.fits
    ql 1 SUPA0111876?.fits
    ql 1 SUPA0146974?.fits
    cd /opt/gen2/data/HSC
    hql 1 HSCA05560*.fits
"""
import os
import glob
import time
import math
import multiprocessing

import numpy as np
from astropy.io import fits

from ginga.rv.plugins.Command import Command, CommandInterpreter
from ginga import AstroImage, cmap
from ginga.gw import Plot, Widgets
from ginga.util import iqcalc, plots, wcs, dp, io_fits
from ginga.misc import Bunch

# add any matplotlib colormaps we have lying around
cmap.add_matplotlib_cmaps(fail_on_import_error=False)

from naoj.hsc.hsc_dr import HyperSuprimeCamDR, hsc_ccd_data

from hsc_mosaic import HSC_Mosaicer

class GView(Command):

    category = 'Utils'

    def __init__(self, fv):
        # superclass defines some variables for us, like logger
        super(GView, self).__init__(fv)

        # get GView preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_GView')
        self.settings.add_defaults(font='Courier', fontsize=12,
                                   plots_in_workspace=False,
                                   history_limit=0,
                                   bee0dir=None,
                                   bee1dir=None,
                                   flatdir=None,
                                   datadir=None)
        self.settings.load(onError='silent')

        self._cmdobj = GViewInterpreter(fv, self)

        fv.add_callback('add-channel', self._cmdobj.add_channel_cb)
        ## fv.add_callback('delete-channel', self.delete_channel)
        ## fv.set_callback('channel-change', self.focus_cb)

        self._plot = None
        self._plot_w = None
        self._plot_idx = 0
        self._plots_in_ws = self.settings.get('plots_in_workspace', False)
        self.w = Bunch.Bunch()


    def build_gui(self, container):

        vbox = Widgets.VBox()

        fr = Widgets.Frame("Output:")

        # plots and output
        if not self._plots_in_ws:
            splitter = Widgets.Splitter(orientation='vertical')

            self.nb = Widgets.TabWidget()
            splitter.add_widget(self.nb)

        fontsize = self.settings.get('fontsize', 12)
        font = self.settings.get('font', 'Courier')
        self.msg_font = self.fv.get_font(font, fontsize)

        vbox.add_widget(Widgets.Label("Output:"))
        tw = Widgets.TextArea(wrap=False, editable=False)
        tw.set_font(self.msg_font)
        tw.set_limit(self.settings.get('history_limit', 0))
        self.hist_w = tw

        if not self._plots_in_ws:
            splitter.add_widget(tw)
            fr.set_widget(splitter)
        else:
            fr.set_widget(tw)

        vbox.add_widget(fr, stretch=1)

        # command line
        vbox.add_widget(Widgets.Label("Type command here:"))
        self.cmd_w = Widgets.TextEntry()
        self.cmd_w.set_font(self.msg_font)
        vbox.add_widget(self.cmd_w, stretch=0)
        self.cmd_w.add_callback('activated', self.exec_cmd_cb)

        hbox = Widgets.HBox()
        btn = Widgets.Button('Detach Plot')
        btn.add_callback('activated', self.detach_plot_cb)
        btn.set_tooltip("Detach current plot; next plot will start a new one")
        hbox.add_widget(btn, stretch=0)
        btn = Widgets.Button('Clear Text')
        btn.add_callback('activated', self.clear_text_cb)
        btn.set_tooltip("Clear the imexam output")
        hbox.add_widget(btn, stretch=0)
        hbox.add_widget(Widgets.Label(''), stretch=1)
        vbox.add_widget(hbox, stretch=0)

        btns = Widgets.HBox()
        btns.set_spacing(4)
        btns.set_border_width(4)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn)
        btns.add_widget(Widgets.Label(''), stretch=1)
        vbox.add_widget(btns)

        self._plot = None
        self._plot_w = None
        self._plot_idx = 0

        container.add_widget(vbox, stretch=1)

    def make_new_figure(self, plot):
        self._plot = plot
        self._plot_idx += 1

        name = "Fig %d" % (self._plot_idx)
        group = 10

        pw = Plot.PlotWidget(plot)

        vbox = Widgets.VBox()
        vbox.add_widget(pw, stretch=1)
        hbox = Widgets.HBox()
        hbox.add_widget(Widgets.Label(''), stretch=1)
        btn = Widgets.Button('Close Plot')
        btn.add_callback('activated', lambda w: self.close_plot(name, vbox))
        hbox.add_widget(btn, stretch=0)
        vbox.add_widget(hbox, stretch=0)

        self._plot_w = vbox

        if self._plots_in_ws:
            ws = self.fv.get_current_workspace()
            tab = self.fv.ds.add_tab(ws.name, vbox, group, name, name,
                                     data=dict(plot=self._plot))
        else:
            self.nb.add_widget(vbox, name)

    def get_current_plot(self):
        return self._plot

    def close_plot(self, name, child):
        if child == self._plot_w:
            self._plot = None
            self._plot_w = None

        if not self._plots_in_ws:
            self.nb.remove(child)

        return True

    def detach_plot_cb(self, w):
        self._plot = None
        self._plot_w = None

    def clear_text_cb(self, w):
        self.hist_w.clear()

    def __str__(self):
        return 'gview'


class GViewInterpreter(CommandInterpreter):

    def __init__(self, fv, plugin):
        super(GViewInterpreter, self).__init__(fv, plugin)

        self.buffers = Bunch.Bunch()

        self.iqcalc = iqcalc.IQCalc(self.logger)

        # Peak finding parameters and selection criteria
        self.radius = 20
        self.settings = {}
        self.max_side = self.settings.get('max_side', 1024)
        self.radius = self.settings.get('radius', 10)
        self.threshold = self.settings.get('threshold', None)
        self.min_fwhm = self.settings.get('min_fwhm', 2.0)
        self.max_fwhm = self.settings.get('max_fwhm', 50.0)
        self.min_ellipse = self.settings.get('min_ellipse', 0.5)
        self.edgew = self.settings.get('edge_width', 0.01)
        self.show_candidates = self.settings.get('show_candidates', False)
        # Report in 0- or 1-based coordinates
        self.pixel_coords_offset = self.settings.get('pixel_coords_offset',
                                                     0.0)

        self.contour_radius = 10

        self.hsc_dr = HyperSuprimeCamDR(logger=self.logger)

        self.sub_bias = True
        # For flat fielding
        self.flat = {}
        self.flat_dir = '.'
        self.flat_filter = None
        self.use_flat = False

    def get_buffer_name(self, bufname):
        try:
            # hack to handle unknown data types
            bufname = int(bufname)
        except ValueError:
            pass
        bufname = str(bufname)
        return bufname

    def add_channel_cb(self, gvshell, channel):
        fi = channel.fitsimage
        bm = fi.get_bindmap()

        # add a new "zview" mode
        bm.add_mode('z', 'zview', mode_type='locked', msg=None)

        # zview had this kind of zooming function
        bm.map_event('zview', (), 'ms_left', 'zoom_out')
        bm.map_event('zview', (), 'ms_right', 'zoom_in')
        bm.map_event('zview', ('ctrl',), 'ms_left', 'zoom_in')
        #bm.map_event('zview', (), 'sc_scroll', 'zoom_origin')

        # borrow some bindings from pan mode
        bm.map_event('zview', (), 'kp_left', 'pan_left')
        bm.map_event('zview', (), 'kp_right', 'pan_right')
        bm.map_event('zview', (), 'kp_up', 'pan_up')
        bm.map_event('zview', (), 'kp_down', 'pan_down')
        bm.map_event('zview', (), 'kp_s', 'pan_zoom_save')
        bm.map_event('zview', (), 'kp_1', 'pan_zoom_set')

        bm.map_event('zview', (), 'kp_p', 'radial-plot')
        bm.map_event('zview', (), 'kp_r', 'radial-plot')
        fi.set_callback('keydown-radial-plot',
                        self.plot_cmd_cb, self.do_radial_plot,
                        "Radial Profile")
        bm.map_event('zview', (), 'kp_e', 'contour-plot')
        fi.set_callback('keydown-contour-plot',
                        self.plot_cmd_cb, self.do_contour_plot,
                        "Contours")
        bm.map_event('zview', (), 'kp_w', 'surface-plot')
        fi.set_callback('keydown-surface-plot',
                        self.plot_cmd_cb, self.do_surface_plot,
                        "Surface")
        bm.map_event('zview', (), 'kp_g', 'gaussians-plot')
        fi.set_callback('keydown-gaussians-plot',
                        self.plot_cmd_cb, self.do_gaussians_plot,
                        "FWHM")

        # bindings customizations
        bd = fi.get_bindings()
        settings = bd.get_settings()

        # ZVIEW has a faster zoom ratio, by default
        settings.set(scroll_zoom_direct_scale=True)

    ##### COMMANDS #####

    def cmd_rd(self, bufname, path, *args):
        """rd bufname path

        Read file from `path` into buffer `bufname`.  If the buffer does
        not exist it will be created.

        If `path` does not begin with a slash it is assumed to be relative
        to the current working directory.
        """
        bufname = self.get_buffer_name(bufname)

        if not path.startswith('/'):
            path = os.path.join(os.getcwd(), path)
        if bufname in self.buffers:
            self.log("Buffer %s is in use. Will discard the previous data" % (
                bufname))
            image = self.buffers[bufname]
        else:
            # new buffer
            image = AstroImage.AstroImage(logger=self.logger)
            self.buffers[bufname] = image

        self.log("Reading file...(%s)" % (path))
        image.load_file(path)
        # TODO: how to know if there is an error
        self.log("File read")

    def cmd_mmbuf(self, bufname, path, wd=42354, ht=42354):
        """mmbuf bufname path

        Create memory-mapped buffer `bufname` by the path `path`.
        If the buffer does not exist it will be created with size
        `ht`x`wd`.

        If `path` does not begin with a slash it is assumed to be relative
        to the current working directory.
        """
        bufname = self.get_buffer_name(bufname)

        if not path.startswith('/'):
            path = os.path.join(os.getcwd(), path)

        hdrpath = path
        dirname, filename = os.path.split(path)
        filepfx, fileext = os.path.splitext(filename)
        datapath = os.path.join(dirname, filepfx + '.mmap')

        if bufname in self.buffers:
            self.log("buffer %s is in use. Will discard the previous data" % (
                bufname))
            image = self.buffers[bufname]
        else:
            # new buffer
            image = AstroImage.AstroImage(logger=self.logger)
            image.set(name=bufname)
            self.buffers[bufname] = image

        if os.path.exists(hdrpath):
            self.logger.info("reading header file: %s" % (hdrpath))
            image.load_file(hdrpath)
            wd, ht = image.get_keywords_list('WIDTH', 'HEIGHT')
        else:
            self.logger.info("writing header file: %s" % (hdrpath))
            wd, ht = int(wd), int(ht)
            hdu = fits.PrimaryHDU()
            header = hdu.header
            header['WIDTH'] = wd
            header['HEIGHT'] = ht
            hdu.writeto(hdrpath)

        if os.path.exists(datapath):
            self.log("reading mmap buffer...(%s)" % (datapath))
            arr_mm = np.memmap(datapath, dtype='float32', mode='w+',
                               shape=(ht, wd))
        else:
            self.log("creating mmap buffer...(%s)" % (datapath))
            arr_mm = np.memmap(datapath, dtype='float32', mode='w+',
                               shape=(ht, wd))
            arr_mm.fill(0.0)

        image.set(statinfo=os.stat(hdrpath), bufdir=dirname, bufpfx=filepfx,
                  hdrpath=hdrpath)

        #image._data = arr_mm
        image.set_data(arr_mm)
        # TODO: how to know if there is an error

        timer = self.fv.get_timer()
        timer.add_callback('expired', self.check_update, image, hdrpath)
        timer.set(3.0)
        self.log("buffer created")

    def check_update(self, timer, image, hdrpath):
        self.logger.debug("checking whether image was modified...")
        old_info = image.get('statinfo', None)
        new_info = os.stat(hdrpath)
        image.set(statinfo=new_info)

        if (old_info is not None) and (new_info.st_mtime > old_info.st_mtime):
            # image is updated
            self.logger.info("image was modified; reloading header...")
            tmp_img = AstroImage.AstroImage(logger=self.logger)
            tmp_img.load_file(hdrpath)

            # replace our header with the new one
            image.set(header=AstroImage.AstroHeader())
            image.update_keywords(tmp_img.get_header())
            wd, ht = image.get_size()
            image.update_keywords(dict(WIDTH=wd, HEIGHT=ht))

            exp_id = image.get_keyword('EXP-ID', None)
            if exp_id is not None:
                self.log("Loaded %s" % (exp_id))

            self.fv.gui_do(image.make_callback, 'modified')

        # check again later
        timer.set(3.0)

    def cmd_tv(self, bufname, *args):
        """tv bufname [min max] [colormap]

        Display buffer `bufname` in the current viewer.  If no viewer
        exists one will be created.

        Optional:
        `min` and `max` specify lo/hi cut levels to scale the image
        data for display.

        `colormap` specifies a color map to use for the image.
        """
        bufname = self.get_buffer_name(bufname)

        if not bufname in self.buffers:
            self.log("!! No such buffer: '%s'" % (bufname))
            return
        image = self.buffers[bufname]

        view = self.make_viewer(bufname)
        view.add_image(image)

        gw = view.viewer

        args = list(args)

        locut = None
        if len(args) > 0:
            try:
                locut = float(args[0])
                hicut = float(args[1])
                args = args[2:]
            except ValueError:
                pass

        if locut is not None:
            gw.cut_levels(locut, hicut)

        if 'bw' in args:
            # replace "bw" with gray colormap
            i = args.index('bw')
            args[i] = 'gray'

        if len(args) > 0:
            cm_name = args[0]
            if cm_name == 'inv':
                gw.invert_cmap()
            else:
                gw.set_color_map(cm_name)

    def cmd_head(self, bufname, *args):
        """head buf [kwd ...]

        List the headers for the image in the named buffer.
        """
        bufname = self.get_buffer_name(bufname)

        if bufname not in self.buffers:
            self.log("No such buffer: '%s'" % (bufname))
            return

        image = self.buffers[bufname]
        header = image.get_header()
        res = []
        # TODO: include the comments
        if len(args) > 0:
            for kwd in args:
                if not kwd in header:
                    res.append("%-8.8s  -- NOT FOUND IN HEADER --" % (kwd))
                else:
                    res.append("%-8.8s  %s" % (kwd, str(header[kwd])))
        else:
            for kwd in header.keys():
                res.append("%-8.8s  %s" % (kwd, str(header[kwd])))

        self.log('\n'.join(res))

    def cmd_exps(self, n=20, hdrs=None):
        """exps  [n=20, time=]

        List the last n exposures in the current directory
        """
        cwd = os.getcwd()
        files = glob.glob(cwd + '/HSCA*[0,2,4,6,8]00.fits')
        files.sort()

        n = int(n)
        files = files[-n:]

        res = []
        for filepath in files:
            with fits.open(filepath, 'readonly', memmap=False) as in_f:
                header = in_f[0].header
            line = "%(EXP-ID)-12.12s  %(HST-STR)12.12s  %(OBJECT)14.14s  %(FILTER01)8.8s" % header

            # add user specified headers
            if hdrs is not None:
                for kwd in hdrs.split(','):
                    fmt = "%%(%s)12.12s" % kwd
                    line += '  ' + (fmt % header)
            res.append(line)

        self.log('\n'.join(res))

    def cmd_lsb(self):
        """lsb

        List the buffers
        """
        names = list(self.buffers.keys())
        names.sort()

        if len(names) == 0:
            self.log("No buffers")
            return

        res = []
        for name in names:
            d = self.get_buffer_info(name)
            d.size = "%dx%d" % (d.width, d.height)
            res.append("%(name)-10.10s  %(size)13s  %(path)s" % d)
        self.log("\n".join(res))

    def cmd_rmb(self, *args):
        """rmb NAME ...

        Remove buffer NAME
        """
        for name in args:
            if name in self.buffers:
                del self.buffers[name]
            else:
                self.log("No such buffer: '%s'" % (name))
        self.cmd_lsb()

    def cmd_rm(self, *args):
        """command to be deprecated--use 'rmb'
        """
        self.log("warning: this command will be deprecated--use 'rmb'")
        self.cmd_rmb(*args)

    def _ql(self, bufname, exp_arg, dr):
        bufname = self.get_buffer_name(bufname)

        if bufname not in self.buffers:
            self.log("!! No such buffer: '%s'" % (bufname))
            return

        image = self.buffers[bufname]

        config = self.plugin.settings
        self.mosaicer = HSC_Mosaicer(self.logger, config)
        self.fv.error_wrap(self.mosaicer.mosaic_exp, exp_arg, image,
                           fork=True)

    def cmd_hql(self, bufname, exp_arg, *args):
        """hql bufname exp_arg
        """
        # arg is an object name or "visit" number
        self._ql(bufname, exp_arg, self.hsc_dr)

    def cmd_bias(self, *args):
        """bias on | off
        """
        if len(args) == 0:
            self.log("bias %s" % (self.sub_bias))
            return
        res = str(args[0]).lower()
        if res in ('y', 'yes', 't', 'true', '1', 'on'):
            self.sub_bias = True
        elif res in ('n', 'no', 'f', 'false', '0', 'off'):
            self.sub_bias = False
        else:
            self.log("Don't understand parameter '%s'" % (onoff))

    def cmd_flat(self, *args):
        """flat on | off
        """
        if len(args) == 0:
            self.log("flat %s" % (self.use_flat))
            return
        res = str(args[0]).lower()
        if res in ('y', 'yes', 't', 'true', '1', 'on'):
            self.use_flat = True
        elif res in ('n', 'no', 'f', 'false', '0', 'off'):
            self.use_flat = False
        else:
            self.log("Don't understand parameter '%s'" % (onoff))

    def cmd_flatdir(self, *args):
        """flatdir /some/path/to/flats
        """
        if len(args) > 0:
            path = str(args[0])
            if not os.path.isdir(path):
                self.log("Not a directory: %s" % (path))
                return
            self.flat_dir = path
        self.log("using (%s) for flats" % (self.flat_dir))

    def get_buffer_info(self, name):
        image = self.buffers[name]
        path = image.get('path', "None")
        res = Bunch.Bunch(dict(name=name, path=path, width=image.width,
                               height=image.height))
        return res

    def make_viewer(self, name):
        if self.fv.has_channel(name):
            channel = self.fv.get_channel(name)
        else:
            channel = self.fv.add_channel(name, num_images=0)

        return channel

    ##### PLOTS #####

    def initialize_plot(self):
        plot = self.plugin.get_current_plot()
        if plot is None:
            plot = plots.Plot(logger=self.logger,
                              width=600, height=600)
            self.plugin.make_new_figure(plot)
        return plot

    def plot_cmd_cb(self, viewer, event, data_x, data_y, fn, title):
        try:
            fn(viewer, event, data_x, data_y)

            #self._plot_w.set_title(title)
            #self._plot_w.raise_()
        except Exception as e:
            self.fv.show_error(e)

        finally:
            # this keeps the focus on the viewer widget, in case a new
            # window was popped up
            #viewer.get_widget().focus()
            pass

    def make_contour_plot(self):
        plot = self.initialize_plot()

        fig = plot.get_figure()
        fig.clf()

        # Replace plot with Contour plot
        plot = plots.ContourPlot(logger=self.logger,
                                 figure=fig,
                                 width=600, height=600)
        kwargs = {'facecolor' if plots.MPL_GE_2_0 else 'axisbg': 'white'}
        plot.add_axis(**kwargs)
        return plot

    def do_contour_plot(self, viewer, event, data_x, data_y):
        self.log("d> (contour plot)", w_time=True)
        try:
            results = self.find_objects(viewer, data_x, data_y)
            qs = results[0]
            x, y = int(round(qs.objx)), int(round(qs.objy))

        except Exception as e:
            self.log("No objects found")
            # we can still proceed with a contour plot at the point
            # where the key was pressed
            x, y = int(round(data_x)), int(round(data_y))

        plot = self.make_contour_plot()
        plot.interpolation = 'nearest'

        image = viewer.get_image()
        data, x1, y1, x2, y2 = image.cutout_radius(x, y,
                                                   self.contour_radius)
        x, y = x - x1, y - y1

        flips = viewer.get_transforms()
        if flips[0]:
            data = np.fliplr(data)
        if flips[1]:
            data = np.flipud(data)

        plot.plot_contours_data(x, y, data, num_contours=12)

        self.fv.ds.raise_tab('GView')
        return True


    def make_gaussians_plot(self):
        plot = self.initialize_plot()

        fig = plot.get_figure()
        fig.clf()

        # Replace plot with FWHM gaussians plot
        plot = plots.FWHMPlot(logger=self.logger,
                              figure=fig,
                              width=600, height=600)
        kwargs = {'facecolor' if plots.MPL_GE_2_0 else 'axisbg': 'white'}
        plot.add_axis(**kwargs)
        return plot

    def do_gaussians_plot(self, viewer, event, data_x, data_y):
        self.log("d> (gaussians plot)", w_time=True)
        try:
            results = self.find_objects(viewer, data_x, data_y)
            qs = results[0]

        except Exception as e:
            self.log("No objects found")
            return

        plot = self.make_gaussians_plot()

        image = viewer.get_image()
        x, y = int(round(qs.objx)), int(round(qs.objy))

        plot.plot_fwhm(x, y, self.radius, image)

        self.fv.ds.raise_tab('GView')
        return True

    def make_radial_plot(self):
        plot = self.initialize_plot()

        fig = plot.get_figure()
        fig.clf()

        # Replace plot with Radial profile plot
        plot = plots.RadialPlot(logger=self.logger,
                                figure=fig,
                                width=700, height=600)
        kwargs = {'facecolor' if plots.MPL_GE_2_0 else 'axisbg': 'white'}
        plot.add_axis(**kwargs)
        return plot

    def do_radial_plot(self, viewer, event, data_x, data_y):
        self.log("d> (radial plot)", w_time=True)
        try:
            results = self.find_objects(viewer, data_x, data_y)
            qs = results[0]
            x, y = int(round(qs.objx)), int(round(qs.objy))

        except Exception as e:
            self.log("No objects found")
            return

        plot = self.make_radial_plot()

        image = viewer.get_image()

        plot.plot_radial(x, y, self.radius, image)

        rpt = self.make_report(image, qs)
        self.log("seeing size %5.2f" % (rpt.starsize))
        # TODO: dump other stats from the report

        # write seeing measurement in upper right corner
        ax = plot.ax
        ax.text(0.75, 0.85, "seeing: %5.2f" % (rpt.starsize),
                bbox=dict(facecolor='green', alpha=0.4, pad=6),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=12)
        plot.draw()

        self.fv.ds.raise_tab('GView')
        return True

    def make_surface_plot(self):
        plot = self.initialize_plot()

        fig = plot.get_figure()
        fig.clf()

        # Replace plot with Surface plot
        plot = plots.SurfacePlot(logger=self.logger,
                                 figure=fig,
                                 width=700, height=600)
        plot.fontsize = 14
        return plot

    def do_surface_plot(self, viewer, event, data_x, data_y):
        self.log("d> (surface plot)", w_time=True)
        try:
            results = self.find_objects(viewer, data_x, data_y)
            qs = results[0]
            x, y = int(round(qs.objx)), int(round(qs.objy))

        except Exception as e:
            self.log("No objects found")
            # we can still proceed with a surface plot at the point
            # where the key was pressed
            x, y = int(round(data_x)), int(round(data_y))

        plot = self.make_surface_plot()

        image = viewer.get_image()

        radius = self.contour_radius * 2
        plot.plot_surface(x, y, radius, image)

        self.fv.ds.raise_tab('GView')
        return True

    def find_objects(self, viewer, x, y):
        #x, y = viewer.get_last_data_xy()
        image = viewer.get_image()

        msg, results, qs = None, [], None
        try:
            data, x1, y1, x2, y2 = image.cutout_radius(int(round(x)),
                                                       int(round(y)),
                                                       self.radius)

            # Find bright peaks in the cutout
            self.logger.info("Finding bright peaks in cutout")
            peaks = self.iqcalc.find_bright_peaks(data,
                                                  threshold=self.threshold,
                                                  radius=self.radius)
            num_peaks = len(peaks)
            if num_peaks == 0:
                raise Exception("Cannot find bright peaks")

            # Evaluate those peaks
            self.logger.info("Evaluating %d bright peaks..." % (num_peaks))
            objlist = self.iqcalc.evaluate_peaks(peaks, data,
                                                 fwhm_radius=self.radius)

            num_candidates = len(objlist)
            if num_candidates == 0:
                raise Exception("Error evaluating bright peaks: no candidates found")

            self.logger.info("Selecting from %d candidates..." % (num_candidates))
            height, width = data.shape
            results = self.iqcalc.objlist_select(objlist, width, height,
                                                 minfwhm=self.min_fwhm,
                                                 maxfwhm=self.max_fwhm,
                                                 minelipse=self.min_ellipse,
                                                 edgew=self.edgew)
            if len(results) == 0:
                raise Exception("No object matches selection criteria")

            # add back in offsets from cutout to result positions
            for qs in results:
                qs.x += x1
                qs.y += y1
                qs.objx += x1
                qs.objy += y1

        except Exception as e:
            msg = str(e)
            self.logger.error("Error finding object: %s" % (msg))
            raise e

        return results

    def make_report(self, image, qs):
        d = Bunch.Bunch()
        try:
            x, y = qs.objx, qs.objy
            equinox = float(image.get_keyword('EQUINOX', 2000.0))

            try:
                ra_deg, dec_deg = image.pixtoradec(x, y, coords='data')
                ra_txt, dec_txt = wcs.deg2fmt(ra_deg, dec_deg, 'str')

            except Exception as e:
                self.logger.warning("Couldn't calculate sky coordinates: %s" % (str(e)))
                ra_deg, dec_deg = 0.0, 0.0
                ra_txt = dec_txt = 'BAD WCS'

            # Calculate star size from pixel pitch
            try:
                header = image.get_header()
                ((xrot, yrot),
                 (cdelt1, cdelt2)) = wcs.get_xy_rotation_and_scale(header)

                starsize = self.iqcalc.starsize(qs.fwhm_x, cdelt1,
                                                qs.fwhm_y, cdelt2)
            except Exception as e:
                self.logger.warning("Couldn't calculate star size: %s" % (str(e)))
                starsize = 0.0

            rpt_x = x + self.pixel_coords_offset
            rpt_y = y + self.pixel_coords_offset

            # make a report in the form of a dictionary
            d.setvals(x = rpt_x, y = rpt_y,
                      ra_deg = ra_deg, dec_deg = dec_deg,
                      ra_txt = ra_txt, dec_txt = dec_txt,
                      equinox = equinox,
                      fwhm = qs.fwhm,
                      fwhm_x = qs.fwhm_x, fwhm_y = qs.fwhm_y,
                      ellipse = qs.elipse, background = qs.background,
                      skylevel = qs.skylevel, brightness = qs.brightness,
                      starsize = starsize,
                      time_local = time.strftime("%Y-%m-%d %H:%M:%S",
                                                 time.localtime()),
                      time_ut = time.strftime("%Y-%m-%d %H:%M:%S",
                                              time.gmtime()),
                      )
        except Exception as e:
            self.logger.error("Error making report: %s" % (str(e)))

        return d

    def __str__(self):
        return 'gview'

#END
