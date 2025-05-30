#
# PFS_Focus.py -- Region selection plugin for fits viewer
#
"""Do a focus fitting for PFS images.

**Plugin Type: Local**

``PFS_Focus`` is a local plugin, which means it is associated with a channel.
An instance can be opened for each channel.

**Usage**

Start the plugin on one of the PFS guide image channels when you are ready
to do a focus fitting.  After a guide image is received in the channel, move
the left and right object selectors (circles) to the objects you which to
measure--one on the left, and one on the right.  Then, click the
"Clear Table" button.

Use the PFS_FOCUSAG observation script to set the Z value using your best
guess for the focus Z value.  After a guide image comes into the channel,
it should be measured and an entry automatically added to the table
(if the values are not automatically added, used the "Measure" button to
add the current measurements to the table).

**Editing the Table**

If you don't want one of the table entries, simply select it and click
"Delete".  If you want to edit it, double click it to copy the values to
the cells below the table.  To remove left or right size measurements,
remove the value from the appropriate "LSize" or "RSize" box--it will be
recorded as "None". When satisfied with the values in the cells, click
"Add Table" to update the values in the table.

Repeat as many times as desired for this Z focus position.  Then, change the
telescope Z value and repeat again.

Cycle through all Z focus positions, populating the table.  Finally, click
"Plot" to do focus fitting on both sides and to calculate the average
focus position.

**The Object Selector**

You can change the object selector size by clicking within the circle and
then dragging the little green dot outward or inward.  The area within the
bounding box of the circle is searched for the center of the object.  If
an object is found, the "X" in the middle of the circle is cyan colored;
if there is a problem determining the center and performing the size
analysis the color will be red.  Selectors can only be dragged within their
respective side of the image.

"""
import math
import time
from collections import OrderedDict
import numpy as np

from ginga.misc import Bunch
from ginga.util import wcs, plots
from ginga.gw import Widgets, Plot
from ginga import GingaPlugin
from ginga.util import iqcalc

from fitsview.util.polynomial import QuadraticFunction

# pixel scale of the PFS guide cameras
pfs_guider_px_scale_X_px_deg = 0.00003944
pfs_guider_px_scale_Y_px_deg = 0.00003694


class PFS_Focus(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(PFS_Focus, self).__init__(fv, fitsimage)

        self.layertag = 'qdas-pfs-focus'

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_PFS_Focus')
        self.settings.add_defaults(sortable=False,
                                   color_alternate_rows=True)
        self.settings.load(onError='silent')

        self.dc = fv.get_draw_classes()
        canvas = self.dc.DrawingCanvas()

        canvas.enable_draw(False)
        canvas.enable_edit(True)
        canvas.set_surface(self.fitsimage)
        canvas.register_for_cursor_drawing(self.fitsimage)
        canvas.set_draw_mode('edit')
        canvas.add_callback('edit-event', self.edit_cb)
        self.canvas = canvas

        lcx, lcy = 266.5, 520
        bbox = self.dc.Box(lcx, lcy, 235, 505, color='yellow')
        bbox.editable = False
        pt = self.dc.Point(lcx, lcy, radius=15, color='red')
        pt.editable = False
        cc = self.dc.Circle(lcx, lcy, radius=30, color='yellow',
                            linewidth=2)
        cc.editable = True
        obj = self.dc.CompoundObject(cc, pt)
        self.left = Bunch.Bunch(cx=lcx, cy=lcy, bbox=bbox, obj=obj)

        rcx, rcy = 804.5, 520
        bbox = self.dc.Box(rcx, rcy, 235, 505, color='blue')
        bbox.editable = False
        pt = self.dc.Point(rcx, rcy, radius=15, color='red')
        pt.editable = False
        cc = self.dc.Circle(rcx, rcy, radius=30, color='blue',
                            linewidth=2)
        cc.editable = True
        obj = self.dc.CompoundObject(cc, pt)
        self.right = Bunch.Bunch(cx=rcx, cy=rcy, bbox=bbox, obj=obj)

        # For image FWHM type calculations
        self.iqcalc = iqcalc.IQCalc(self.logger)

        self.radius = 10
        self.threshold = None
        self.use_new_algorithm = True

        self.columns = [('AG Image', 'imname'),
                        ('Z', 'z'),
                        ('L Size', 'lsize'),
                        ('R Size', 'rsize'),
                        ('L-R Diff', 'lmrdiff'),
                        ]
        self.tree_dict = OrderedDict()

        self.qf_l = QuadraticFunction(logger=self.logger)
        self.qf_r = QuadraticFunction(logger=self.logger)

        self.gui_up = False

    def build_gui(self, container, future=None):

        top = Widgets.VBox()
        top.set_border_width(0)
        orientation = 'vertical'

        paned = Widgets.Splitter(orientation=orientation)
        self.w.splitter = paned

        vbox = Widgets.VBox()
        vbox.set_border_width(2)

        table = Widgets.TreeView(auto_expand=True,
                                 selection='single',
                                 use_alt_row_color=True)
        self.w.table = table
        table.setup_table(self.columns, 1, 'imname')
        table.add_callback('activated', self.edit_table_cb)

        vbox.add_widget(table, stretch=1)

        captions = (("ImName:", 'label', 'imname', 'entry',
                     "Z:", 'label', 'Z', 'entry',
                     "LSize:", 'label', 'LSize', 'entry',
                     "RSize:", 'label', 'RSize', 'entry'),
                    ('Add Table', 'button',
                     'Delete', 'button', 'Clear Table', 'button',
                     'Auto', 'checkbox', 'Measure', 'button'),
                    ('Plot', 'button',
                     ### 'Debug', 'combobox', # Comment-out for now
                     ),
                    )

        w, b = Widgets.build_info(captions, orientation=orientation)
        self.w.update(b)
        b.measure.add_callback('activated',
                               lambda w: self.update_values())
        b.measure.set_tooltip("Measure current image and add to table")
        b.add_table.add_callback('activated',
                                 lambda w: self.add_to_table())
        b.add_table.set_tooltip("Add/edit updated cells to table")
        b.delete.add_callback('activated',
                              lambda w: self.delete_entries())
        b.delete.set_tooltip("Delete selected entry from table")
        b.clear_table.add_callback('activated',
                                   lambda w: self.clear_table())
        b.clear_table.set_tooltip("Clear the entire table")
        b.plot.add_callback('activated',
                            lambda w: self.plot_data())
        b.plot.set_tooltip("Plot the focus fitting curves from table")
        b.auto.set_tooltip("Automatically add measurements to table")

        # for i in range(0, 8):
        #     b.debug.append_text(str(i))
        # b.debug.add_callback('activated', self._debug_set_table)

        vbox.add_widget(w, stretch=0)
        paned.add_widget(vbox)

        self.focus_plot = plots.Plot(logger=self.logger,
                                     width=400, height=400)
        self.focus_plot.add_axis(facecolor='white')
        self.init_plot()

        self.lmr_plot = plots.Plot(logger=self.logger,
                                   width=400, height=400)
        self.lmr_plot.add_axis(facecolor='white')
        self.init_lmr_plot()

        self.w.plot_w = Plot.PlotWidget(self.focus_plot)
        self.w.plot_w.resize(400, 400)

        self.w.lmr_plot_w = Plot.PlotWidget(self.lmr_plot)
        self.w.lmr_plot_w.resize(400, 400)

        self.tabWidget = Widgets.TabWidget()
        self.tabWidget.add_widget(self.w.plot_w, title='Focus')
        self.tabWidget.add_widget(self.w.lmr_plot_w, title='L-R')

        paned.add_widget(self.tabWidget)

        top.add_widget(paned, stretch=1)
        #vbox.add_widget(Widgets.Label(''), stretch=1)

        captions = (("Best Focus:", 'label', 'best_z', 'llabel'),
                    )

        w, b = Widgets.build_info(captions, orientation=orientation)
        self.w.update(b)
        top.add_widget(w, stretch=0)
        b.best_z.set_text('')

        btns = Widgets.HBox()
        btns.set_border_width(4)
        btns.set_spacing(4)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn)
        #btn = Widgets.Button("Help")
        #btn.add_callback('activated', lambda w: self.help())
        #btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)

        container.add_widget(top, stretch=1)

        self.gui_up = True

    def withdraw_qdas_layers(self):
        p_canvas = self.fitsimage.get_canvas()
        tags = p_canvas.get_tags_by_tag_pfx('qdas-')
        for tag in tags:
            try:
                p_canvas.delete_object_by_tag(tag)
            except KeyError:
                pass

    def start(self):
        # remove all qdas canvases
        self.withdraw_qdas_layers()

        # insert our canvas to fitsimage if it is not already
        p_canvas = self.fitsimage.get_canvas()
        try:
            obj = p_canvas.get_object_by_tag(self.layertag)

        except KeyError:
            # Add canvas layer
            p_canvas.add(self.canvas, tag=self.layertag)

        self.canvas.delete_all_objects(redraw=False)

        # draw our two bounding boxes
        self.canvas.add(self.left.bbox, redraw=False)
        self.canvas.add(self.right.bbox, redraw=False)
        # add detected objects
        self.canvas.add(self.left.obj, redraw=False)
        self.canvas.add(self.right.obj, redraw=False)
        self.fitsimage.redraw(whence=3)

        self.resume()

    def pause(self):
        self.logger.debug("disabling canvas")
        self.canvas.ui_set_active(False)

    def resume(self):
        # turn off any mode user may be in
        self.modes_off()

        self.canvas.ui_set_active(True)

    def stop(self):
        self.pause()
        # remove our canvas from image
        p_canvas = self.fitsimage.get_canvas()
        try:
            p_canvas.delete_object_by_tag(self.layertag)

        except KeyError:
            pass
        self.gui_up = False

    def close(self):
        chname = self.fv.get_channel_name(self.fitsimage)
        self.fv.stop_local_plugin(chname, str(self))
        return True

    def get_info(self):
        # TODO: this will eventually fetch the current Z value from Gen2
        stat_d = {'TSCL.Z': None}
        self.logger.info("calling fetch!")
        try:
            self.fv.call_global_plugin_method('Gen2Int', 'fetch', [stat_d], {})
            self.logger.info("fetch returned with %s" % str(stat_d))
            z = float(stat_d['TSCL.Z'])
        except Exception as e:
            self.logger.error("failed calling status: {}".format(e),
                              exc_info=True)
            z = None
        return Bunch.Bunch(z=z)

    def update_values(self):
        self._redo()

    def add_to_table(self):
        imname = self.w.imname.get_text().strip()
        lsize = self.w.lsize.get_text().strip()
        if len(lsize) == 0 or lsize.lower() == 'none':
            lsize = None
        else:
            lsize = float(lsize)
        rsize = self.w.rsize.get_text().strip()
        if len(rsize) == 0 or rsize.lower() == 'none':
            rsize = None
        else:
            rsize = float(rsize)
        lmrdiff = lsize - rsize if lsize is not None and rsize is not None else None

        z = self.w.z.get_text().strip()
        if len(z) == 0 or z.lower() == 'none':
            z = None
        else:
            z = float(z)
        bnch = Bunch.Bunch(imname=imname, z=z,
                           lsize=lsize, rsize=rsize, lmrdiff=lmrdiff)

        # add to table
        self.tree_dict[imname] = bnch
        self.w.table.set_tree(self.tree_dict)
        self.w.table.scroll_to_end()
        #self.w.table.set_optimal_column_widths()

    def edit_table_cb(self, widget, selection):
        bnch = list(selection.values())[0]
        self.w.imname.set_text(bnch.imname)
        self.w.z.set_text("{z:.4f}".format(z=bnch.z))
        self.w.lsize.set_text("{lsize:.4f}".format(lsize=bnch.lsize))
        self.w.rsize.set_text("{rsize:.4f}".format(rsize=bnch.rsize))

    def delete_entries(self):
        """Delete the selected entries from the table."""
        selected = self.w.table.get_selected()
        tree_dict = OrderedDict()
        for key, val in self.tree_dict.items():
            if key not in selected:
                tree_dict[key] = val

        self.tree_dict = tree_dict
        self.w.table.set_tree(self.tree_dict)

    def clear_table(self):
        """Clear table and plot."""
        self.init_plot()

        self.tree_dict = OrderedDict()
        self.w.table.set_tree(self.tree_dict)
        self.w.best_z.set_text('')

    def init_plot(self):
        # set up matplotlib axis here as needed
        ax = self.focus_plot.get_axis()
        ax.cla()
        ax.grid(True)
        ax.set_title("PFS Focus Fitting",
                     fontdict=dict(fontsize=14, fontweight='bold',
                                   color='darkgreen'))
        ax.set_xlabel("Z Position",
                      fontdict=dict(fontsize=12, fontweight='bold'))
        ax.set_ylabel("Star Size",
                      fontdict=dict(fontsize=12, fontweight='bold'))

        self.focus_plot.draw()

    def init_lmr_plot(self):
        # set up matplotlib axis here as needed
        ax = self.lmr_plot.get_axis()
        ax.cla()
        ax.grid(True)
        ax.set_title("PFS Focus L - R Diff",
                     fontdict=dict(fontsize=14, fontweight='bold',
                                   color='darkgreen'))
        ax.set_xlabel("Z Position",
                      fontdict=dict(fontsize=12, fontweight='bold'))
        ax.set_ylabel("L - R",
                      fontdict=dict(fontsize=12, fontweight='bold'))

        self.lmr_plot.draw()

    def plot_data(self):
        self.plot_curves()
        self.plot_lmr()

    def plot_curves(self):
        self.init_plot()
        values = [(val.z, val.lsize, val.rsize)
                  for val in self.tree_dict.values()]

        values = np.array(values, dtype=float)
        x_points, lsize, rsize = values.T
        #print(f'x={x_points}, l={lsize}, r={rsize}')
        self.logger.info(f'x_points {x_points}')
        self.logger.info(f'lsize {lsize}')
        self.logger.info(f'rsize {rsize}')

        lidx = np.isfinite(x_points) & np.isfinite(lsize)
        ridx = np.isfinite(x_points) & np.isfinite(rsize)

        # figure
        fig = self.focus_plot.get_figure()
        ax = self.focus_plot.get_axis()
        # plot curves

        ax.plot(x_points, lsize, 'yo')
        ax.plot(x_points, rsize, 'bo')

        def vertex(x, y):
            bbox_args = dict(boxstyle="round", fc="cyan", alpha=0.1)
            ax.annotate('Vertex(%.2f, %.2f)'%(x,y)  , xy=(x,y), xytext=(x, y),
                               size=20, color='g',
                               bbox=bbox_args,
                               ha='center',
                               va='top')

        def err_msg(msg, x, y):
            bbox_args = dict(boxstyle="round", fc="red", alpha=0.6)
            ax.annotate(msg, xy=(x, y), xytext=(x, y),
                        size=20, color='g',
                        bbox=bbox_args,
                        ha='center',
                        va='bottom')

        errors = []
        left = "LEFT"
        right = "RIGHT"

        xs = len(x_points) * 10 # 10 times more x points
        x_more_points = np.linspace(x_points[0], x_points[-1], xs)

        try:
            self.qf_l.coefficient(x_points[lidx], lsize[lidx])
            self.logger.info(f'qf_l coeffs={self.qf_l.coeffs}')
            qf_l_func = self.qf_l.quadratic()
            ax.plot(x_more_points, qf_l_func(x_more_points), 'y-', linewidth=2)

            l_min_x, l_min_y = self.qf_l.min_vertex()
            #print(f'lside. l_min_x={l_min_x}, l_min_y={l_min_y}')
        except Exception as e:
            #print(f'qf lsize to do {e}')
            errors.append(left)

        try:
            self.qf_r.coefficient(x_points[ridx], rsize[ridx])
            self.logger.info(f'qf_r coeffs={self.qf_r.coeffs}')
            qf_r_func = self.qf_r.quadratic()
            ax.plot(x_more_points, qf_r_func(x_more_points), 'b-', linewidth=2)

            r_min_x, r_min_y = self.qf_r.min_vertex()
            #print(f'rside. r_min_x={r_min_x}, r_min_y={r_min_y}')
        except Exception as e:
            errors.append(right)


        if len(errors) == 2:
            x = x_points.mean()
            y = 0
            msg = "L/R coefficient calc errors"
            err_msg(msg, x, y)
            self.logger.error(f'{msg}')

        elif any(left in e for e in errors):
            x = r_min_x
            y = r_min_y
            vertex(x, y)

            msg = "Left coefficient calc error"
            y = rsize[ridx][-1]
            err_msg(msg, x, y)
            self.logger.error(f'{msg}')

        elif any(right in e for e in errors):
            x = l_min_x
            y = l_min_y
            vertex(x, y)

            msg = "Right coefficient calc error"
            y = lsize[lidx][-1]
            err_msg(msg, x, y)
            self.logger.error(f'{msg}')

        else:
            # determine best Z position as mean of left and right minimums
            y = (l_min_y + r_min_y) * 0.5
            x = (l_min_x + r_min_x) * 0.5

            vertex(x, y)
            self.logger.info(f'l_min_x {l_min_x} l_min_y {l_min_y}')
            self.logger.info(f'r_min_x {r_min_x} r_min_y {r_min_y}')
            self.logger.info(f'x {x} y {y}')

        best_z = x
        self.w.best_z.set_text(f"Z = {best_z:.4f}")
        self.logger.info(f'best_z {best_z}')

        self.focus_plot.draw()

    def plot_lmr(self):
        self.init_lmr_plot()
        values = [(val.z, val.lmrdiff)
                  for val in self.tree_dict.values()]

        values = np.array(values, dtype=float)
        x_points, lmrdiff = values.T
        #print(f'x={x_points}, lmr={lmrdiff}')
        self.logger.info(f'x_points {x_points}')
        self.logger.info(f'lmrdiff {lmrdiff}')

        lmridx = np.isfinite(x_points) & np.isfinite(lmrdiff)

        # figure
        fig = self.lmr_plot.get_figure()
        ax = self.lmr_plot.get_axis()
        # plot curves

        ax.plot(x_points, lmrdiff, 'yo')

        def zero_cross(zcross, x, y):
            bbox_args = dict(boxstyle="round", fc="cyan", alpha=0.1)
            ax.annotate('Zero Cross @ Z=%.2f'%(zcross), xy=(x,y), xytext=(x, y),
                               size=20, color='g',
                               bbox=bbox_args,
                               ha='center',
                               va='top')

        def err_msg(msg, x, y):
            bbox_args = dict(boxstyle="round", fc="red", alpha=0.6)
            ax.annotate(msg, xy=(x, y), xytext=(x, y),
                        size=20, color='g',
                        bbox=bbox_args,
                        ha='center',
                        va='bottom')

        xs = len(x_points) * 10 # 10 times more x points
        x_more_points = np.linspace(x_points[0], x_points[-1], xs)

        try:
            coeffs = np.polyfit(x_points[lmridx], lmrdiff[lmridx], 1)
            fit_fn = np.poly1d(coeffs)
            self.logger.info(f'coeffs {coeffs}')
            if math.isclose(coeffs[0], 0.0, abs_tol=0.001):
                msg = 'Slope close to 0'
                x = x_points.mean()
                y = lmrdiff[lmridx][-1]
                err_msg(msg, x, y)
                self.logger.error(f'{msg}')
            else:
                zcross = -coeffs[1] / coeffs[0]
                x = x_points.mean()
                y = lmrdiff[lmridx][-1]
                self.logger.info(f'zcross {zcross} y {y}')
                zero_cross(zcross, x, y)
            ax.plot(x_more_points, fit_fn(x_more_points), 'y-', linewidth=2)
        except Exception as e:
            msg = 'Error computing linear fit'
            x = x_points.mean()
            y = 0
            err_msg(msg, x, y)
            self.logger.error(f'{msg}: {str(e)}')

        self.lmr_plot.draw()

    def _redo(self):
        # get image name
        #image = self.fitsimage.get_image()
        #imname = image.get('name')
        imname = time.strftime("AG%H:%M:%S", time.localtime())

        # get status info
        try:
            info = self.get_info()
            z = float(info.z)
        except Exception as e:
            z = None

        # measurements on left and right sides
        p = self.setup_side(self.left.obj.objects[0], self.left.obj,
                            self.left.bbox)
        # for testing/debugging
        #p.starsize += 2 ** np.random.random() * np.abs(z)

        q = self.setup_side(self.right.obj.objects[0], self.right.obj,
                            self.right.bbox)
        # for testing/debugging
        #q.starsize += 2 ** np.random.random() * np.abs(z)

        lsize = p.starsize
        rsize = q.starsize
        lmrdiff = lsize - rsize if lsize is not None and rsize is not None else None
        bnch = Bunch.Bunch(imname=imname, z=z,
                           lsize=p.starsize, rsize=q.starsize, lmrdiff=lmrdiff)

        # add to table
        self.tree_dict[imname] = bnch
        self.w.table.set_tree(self.tree_dict)
        self.w.table.scroll_to_end()

    def redo(self):
        if not self.gui_up:
            return
        if not self.w.auto.get_state():
            return
        self._redo()

    def setup_side(self, bbox_search, obj, bbox_limits):
        x1, y1, x2, y2 = bbox_search.get_llur()
        p = Bunch.Bunch()
        circle = obj.objects[0]
        point = obj.objects[1]
        try:
            image = self.fitsimage.get_image()
            if image is None:
                raise ValueError("No image is present in viewer")

            qualsize = self.iqcalc.qualsize
            qs = qualsize(image, x1, y1, x2, y2,
                          radius=self.radius, threshold=self.threshold,
                          minelipse=0.2)

            xl1, yl1, xl2, yl2 = bbox_limits.get_llur()
            if not (xl1 <= qs.objx <= xl2) or not (yl1 <= qs.objy <= yl2):
                raise ValueError("object center outside of search area")

            point.color = 'cyan'
            # Calculate X/Y of center of star
            circle.x = point.x = qs.objx
            circle.y = point.y = qs.objy

            header = image.get_header()
            rot, cdelt1, cdelt2 = wcs.get_rotation_and_scale(header)
            pscale_x, pscale_y = np.abs(cdelt1), np.abs(cdelt2)
            self.logger.info(f"fwhm_x: {qs.fwhm_x}, fwhm_y: {qs.fwhm_y}")
            self.logger.info(f"px scale measured from WCS ({pscale_x},{pscale_y})")
            starsize = self.iqcalc.starsize(qs.fwhm_x, pscale_x,
                                            qs.fwhm_y, pscale_y)
            self.logger.info(f"WCS computed star size ({starsize})")
            # TEMP: pixel scale encoded in PFS AG images is not correct
            # Manually override with pixel scale provided by Kawanomoto-san
            pscale_x = pfs_guider_px_scale_X_px_deg
            pscale_y = pfs_guider_px_scale_Y_px_deg
            self.logger.info(f"px scale override ({pscale_x},{pscale_y})")
            starsize = self.iqcalc.starsize(qs.fwhm_x, pscale_x,
                                            qs.fwhm_y, pscale_y)
            self.logger.info(f"override computed star size ({starsize})")

            p.fwhm = qs.fwhm
            p.starsize = starsize
            p.brightness = qs.brightness
            p.skylevel = qs.skylevel
            p.obj_x = qs.objx
            p.obj_y = qs.objy

        except Exception as e:
            self.logger.error("Error calculating quality metrics: %s" % (
                str(e)))

            # set region
            x, y = (x1 + x2) * 0.5, (y1 + y2) * 0.5
            circle.x = point.x = x
            circle.y = point.y = y
            point.color = 'red'
            p.starsize = None

        self.canvas.redraw(whence=3)
        return p

    def edit_cb(self, canvas, obj):
        # <-- user moved or edited a circle
        # if the object position is outside our bbox, reset it to the center
        # of the bbox. Otherwise, move the point to the center of the circle
        x, y = obj.x, obj.y
        if obj is self.left.obj.objects[0]:
            x1, y1, x2, y2 = self.left.bbox.get_llur()
            if not (x1 <= x <= x2) or not (y1 <= y <= y2):
                x, y = self.left.bbox.x, self.left.bbox.y
                obj.x = x
                obj.y = y
            self.left.obj.objects[1].x = obj.x
            self.left.obj.objects[1].y = obj.y

            p = self.setup_side(obj, self.left.obj, self.left.bbox)
            if p.starsize is None:
                self.w.lsize.set_text("Error")
            else:
                self.w.lsize.set_text("{lsize:.4f}".format(lsize=p.starsize))

        elif obj is self.right.obj.objects[0]:
            x1, y1, x2, y2 = self.right.bbox.get_llur()
            if not (x1 <= x <= x2) or not (y1 <= y <= y2):
                x, y = self.right.bbox.x, self.right.bbox.y
                obj.x = x
                obj.y = y
            self.right.obj.objects[1].x = obj.x
            self.right.obj.objects[1].y = obj.y

            q = self.setup_side(obj, self.right.obj, self.right.bbox)
            if q.starsize is None:
                self.w.rsize.set_text("Error")
            else:
                self.w.rsize.set_text("{rsize:.4f}".format(rsize=q.starsize))

        #self.redo()

        self.fitsimage.redraw(whence=3)
        return True

    def _debug_set_table(self, w, idx):
        if idx == 1:
            # single data with None per X point
            values = [(-3.0, 4.0, 3.0), (-2.0, None, 2.0), (-1.0, 2.0, 1.0), (0.0, 1.0, 0.0), (1.0, 2.0, None), (2.0, 3.0, 2.0), (3.0, 4.0, 3.0)]

        elif idx == 2:
            # multiple data per X point
            values = [(-3.0, 4.0, 3.0), (-3.0, 4.1, 2.9), (-2.0, 3.0, 2.0), (-2.0, 3.1, 1.9),  (-1.0, 2.0, 1.0),  (-1.0, 2.1, 0.9), (0.0, 1.5, 0.5), (0.0, 1.7, 0.3), (1.0, 2.0, 1.0), (1.0, 2.2, 1.2), (2.0, 3.0, 2.0), (2.0, 2.7, 2.3), (3.0, 4.0, 3.0), (3.0, 4.2, 3.1)]

        elif idx == 3:
            # multiple data with None per X point
            values = [(-3.0, 4.0, 3.0), (-3.0, None, 2.9), (-2.0, 3.0, 2.0), (-2.0, 3.1, 1.9),  (-1.0, 2.0, 1.0),  (-1.0, 2.1, 0.9), (0.0, 1.5, 0.5), (0.0, 1.7, 0.3), (1.0, 2.0, None), (1.0, 2.2, 1.2), (2.0, None, 2.0), (2.0, None, 2.3), (3.0, 4.0, 3.0), (3.0, 4.2, 3.1)]

        elif idx == 4:
            # left side error
            values = [(-3.0, None, 3.0), (-3.0, None, 2.9), (-2.0, None, 2.0), (-2.0, None, 1.9),  (-1.0, None, 1.0),  (-1.0, None, 0.9), (0.0, None, 0.5), (0.0, None, 0.3), (1.0, None, None), (1.0, None, 1.2), (2.0, None, 2.0), (2.0, None, 2.3), (3.0, None, 3.0), (3.0, None, 3.1)]

        elif idx == 5:
            # right side error
            values = [(-3.0, 3.0, None), (-3.0, 2.9, None), (-2.0, 2.0, None), (-2.0, 1.9, None),  (-1.0, 1.0, None),  (-1.0, 0.9, None), (0.0, 0.5, None), (0.0, 0.3, None), (1.0, None, None), (1.0, 1.2, None), (2.0, 2.0, None), (2.0, 2.3, None), (3.0, 3.0, None), (3.0, 3.1, None)]

        elif idx == 6:
            # left/right side errors
            values = [(-3.0, None, None), (-3.0, None, None), (-2.0, None, None), (-2.0, None, None),  (-1.0, None, None),  (-1.0, None, None), (0.0, None, None), (0.0, None, None), (1.0, None, None), (1.0, None, None), (2.0, None, None), (2.0, None, None), (3.0, None, None), (3.0, None, None)]

        elif idx == 7:
            # Real data from HST 23 July 2023 23:24:41
            values = [(2.95, 2.63596, 2.48303), (2.95, 2.5551, 3.041967), (2.95, 2.24997, 2.39144978), (3.05, 2.474475, 2.483757), (3.05, 1.4191977, 1.72046), (3.05, 1.54348, 1.9384089), (3.15, 2.31308, 1.733852796)]

        else:
            # single data per X point
            values = [(-3.0, 4.0, 3.0), (-2.0, 3.0, 2.0), (-1.0, 2.0, 1.0), (0.0, 1.5, 0.5), (1.0, 2.0, 1.0), (2.0, 3.0, 2.0), (3.0, 4.0, 3.0)]

        # DEBUG: set table to match values
        self.clear_table()
        for i, tup in enumerate(values):
            z, lsize, rsize = tup
            lmrdiff = lsize - rsize if lsize is not None and rsize is not None else None
            imname = f"AG{i}"
            self.tree_dict[imname] = Bunch.Bunch(imname=imname, z=z,
                                                 lsize=lsize, rsize=rsize, lmrdiff=lmrdiff)
        self.w.table.set_tree(self.tree_dict)

    def __str__(self):
        return 'pfs_focus'
