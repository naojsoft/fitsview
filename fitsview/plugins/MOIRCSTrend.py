# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
"""
A plugin for generating a plot of the values along a line or path.

**Plugin Type: Global**

``MOIRCSTrend`` is a global plugin, which means only one instance can
be opened.

**Usage**

"""
from types import SimpleNamespace
import dateutil.parser
from dateutil import tz
import numpy as np
import threading

from ginga.gw import Widgets
from ginga import GingaPlugin, colors

try:
    from ginga.gw import Plot
    from ginga.util import plots

    import matplotlib.dates as mpl_dt
    import matplotlib as mpl
    from matplotlib.ticker import FormatStrFormatter

    have_mpl = True
except ImportError:
    have_mpl = False

from g2base.astro.frame import Frame


__all__ = ['MOIRCSTrend']


class MOIRCSTrend(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        super().__init__(fv)

        self._split_sizes = [400, 500]
        self._data = dict(det1=SimpleNamespace(color='green', points={}),
                          det2=SimpleNamespace(color='purple', points={}))
        self.tz = tz.gettz('US/Hawaii')

        # get MOIRCSTrend preferences
        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_MOIRCSTrend')
        self.settings.add_defaults(orientation='vertical',
                                   show_legend=True)
        self.settings.load(onError='silent')

        self.gui_up = False

    def build_gui(self, container):
        if not have_mpl:
            raise ImportError('Install matplotlib to use this plugin')

        top = Widgets.VBox()
        top.set_border_width(4)

        # Make the cuts plot
        box, sw, orientation = Widgets.get_oriented_box(container,
                                                        orientation=self.settings.get('orientation', None))
        box.set_margins(4, 4, 4, 4)
        box.set_spacing(2)

        paned = Widgets.Splitter(orientation=orientation)
        self.w.splitter = paned

        # Add Tab Widget
        nb = Widgets.TabWidget(tabpos='top')
        paned.add_widget(Widgets.hadjust(nb, orientation))

        self.trend_plot = plots.Plot(logger=self.logger,
                                     width=400, height=400)
        self.plot = Plot.PlotWidget(self.trend_plot)
        self.plot.resize(400, 400)
        ax = self.trend_plot.add_axis()
        ax.grid(True)

        paned.add_widget(sw)
        paned.set_sizes(self._split_sizes)

        top.add_widget(paned, stretch=5)

        # Add plot to its tab
        vbox_trend = Widgets.VBox()
        vbox_trend.add_widget(self.plot, stretch=1)
        nb.add_widget(vbox_trend, title="Trend")

        btns = Widgets.HBox()
        btns.set_border_width(4)
        btns.set_spacing(3)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        # btn = Widgets.Button("Help")
        # btn.add_callback('activated', lambda w: self.help())
        # btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)

        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)

        self.gui_up = True

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def start(self):
        self.trend_plot.set_titles(rtitle="MOIRCS Trend")

        if self.gui_up:
            self.fv.gui_do(self.replot_all)
        self.resume()

    def pause(self):
        pass

    def resume(self):
        pass

    def stop(self):
        self.gui_up = False
        self._split_sizes = self.w.splitter.get_sizes()
        # remove the canvas from the image
        self.fv.show_status("")

    def redo(self, channel, image):
        """This is called when a new image arrives or the data in the
        existing image changes.
        """
        path = image.get('path', None)
        if path is None:
            return
        fr = Frame(path)
        if fr.inscode != 'MCS' or fr.frametype != 'A':
            # only MOIRCS A frames
            return
        hdr = image.get_header()
        if 'DET-ID' not in hdr:
            return
        # skip non-image files
        if 'OBS-MOD' not in hdr:
            return
        if hdr['OBS-MOD'].strip() != 'IMAG':
            return
        det_id = int(hdr['DET-ID'])
        frame_no = fr.number
        # get configuration for this detector
        my_cfg = self._data[f'det{det_id}']

        if frame_no in my_cfg.points:
            # we have seen this frame before--don't process it again
            return

        # collect info
        data = image.get_data()
        med = np.nanmedian(data)

        # get obs date/time
        obs_date = hdr['DATE-OBS'].strip()
        obs_time = hdr['UT'].strip()
        dt = dateutil.parser.parse("{} {}".format(obs_date, obs_time)).replace(tzinfo=tz.UTC)
        # add to history
        my_cfg.points[frame_no] = SimpleNamespace(median=med,
                                                  time=dt,
                                                  number=frame_no)
        if self.gui_up:
            self.fv.gui_do(self.replot_all)

    def add_legend(self):
        """Add or update Cuts plot legend."""
        self.trend_plot.ax.legend(['det 1', 'det 2'], loc='best',
                                  shadow=True, fancybox=True,
                                  prop={'size': 8}, labelspacing=0.2)

    def replot_all(self):
        if not self.gui_up:
            return
        det1 = self._data['det1']
        _keys = sorted(det1.points.keys())
        det1_pts = None
        if len(_keys) > 0:
            det1_pts = np.array([(det1.points[k].time.astimezone(self.tz),
                                  det1.points[k].median)
                                 for k in _keys])

        det2 = self._data['det2']
        _keys = sorted(det2.points.keys())
        det2_pts = None
        if len(_keys) > 0:
            det2_pts = np.array([(det2.points[k].time.astimezone(self.tz),
                                  det2.points[k].median)
                                 for k in _keys])

        ax = self.trend_plot.ax
        ax.cla()
        ax.grid()

        # set major ticks to hours
        majorTick = mpl_dt.HourLocator(tz=self.tz)
        majorFmt = mpl_dt.DateFormatter('%Hh', tz=self.tz)
        # set minor ticks to 15 min intervals
        minorTick = mpl_dt.MinuteLocator(list(range(0,59,15)), tz=self.tz)

        if det1_pts is not None:
            x_arr, y_arr = det1_pts.T
            ax.plot_date(x_arr, y_arr, '-', color=self._data['det1'].color,
                         aa=True, tz=self.tz)
        if det2_pts is not None:
            x_arr, y_arr = det2_pts.T
            ax.plot_date(x_arr, y_arr, '-', color=self._data['det2'].color,
                         aa=True, tz=self.tz)

        self.trend_plot.set_titles(title="MOIRCS Trend",
                                   xtitle="Time (HST)", ytitle="Median")
        if self.settings.get('show_legend', False):
            self.add_legend()

        self.trend_plot.draw()

    def __str__(self):
        return 'moircstrend'
