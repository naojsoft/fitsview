#
# AgAreaSelection.py -- Ag area selection plugin for fits viewer
# 
# Eric Jeschke (eric@naoj.org)
#
import gtk
from ginga.gtkw import GtkHelp

from ginga.misc import Bunch

from ginga.gtkw import ImageViewCanvasGtk
from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
from ginga import GingaPlugin

# Local application imports
from util import g2calc

class QDASPlugin(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        super(QDASPlugin, self).__init__(fv, fitsimage)

        self.layertag = 'qdas-agareaselection'

        canvas = CanvasTypes.DrawingCanvas()
        canvas.enable_draw(True)
        self.canvas = canvas

        canvas.set_callback('cursor-down', self.drag)
        canvas.set_callback('cursor-move', self.drag)
        canvas.set_callback('cursor-up', self.update)
        canvas.set_callback('draw-event', self.setpickregion)
        canvas.set_drawtype('rectangle', color='cyan', linestyle='dash',
                            drawdims=True)
        canvas.setSurface(self.fitsimage)

    def withdraw_qdas_layers(self):
        tags = self.fitsimage.getTagsByTagpfx('qdas-')
        for tag in tags:
            try:
                self.fitsimage.deleteObjectByTag(tag)
            except:
                pass
        


class AgAreaSelection(QDASPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(AgAreaSelection, self).__init__(fv, fitsimage)

        self.exptime = 6000
        
        self.pick_qs = None
        self.picktag = None
        self.pickcolor = 'green'
        self.eregcolor = 'pink'

        self.dx = 30
        self.dy = 30
        # this is the maximum size a side can be
        self.max_len = 512
        self.agarea = 'AG1'

        self.iqcalc = g2calc.IQCalc(self.logger)
        self.radius = 10
        self.threshold = None
        self.use_new_algorithm = False
        
    def build_gui(self, container, future=None):
        sw = gtk.ScrolledWindow()
        sw.set_border_width(2)
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

        vbox = gtk.VBox()
        sw.add_with_viewport(vbox)

        self.msgFont = self.fv.getFont("sansFont", 14)
        tw = gtk.TextView()
        tw.set_wrap_mode(gtk.WRAP_WORD)
        tw.set_left_margin(4)
        tw.set_right_margin(4)
        tw.set_editable(False)
        tw.set_left_margin(4)
        tw.set_right_margin(4)
        tw.modify_font(self.msgFont)
        self.tw = tw

        fr = gtk.Frame(" Instructions ")
        fr.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
        fr.set_label_align(0.1, 0.5)
        fr.add(tw)
        vbox.pack_start(fr, padding=4, fill=True, expand=False)
        
        nb = gtk.Notebook()
        #nb.set_group_id(group)
        #nb.connect("create-window", self.detach_page, group)
        nb.set_tab_pos(gtk.POS_BOTTOM)
        nb.set_scrollable(True)
        nb.set_show_tabs(True)
        nb.set_show_border(False)
        self.w.nb2 = nb
        vbox.pack_start(nb, padding=4, fill=True, expand=True)

        captions = (
            ('Object_X', 'label', 'Object_Y', 'label'),
            ('RA', 'label', 'DEC', 'label'), ('Equinox', 'label'),
            ('Sky Level', 'label', 'Brightness', 'label'), 
            ('FWHM', 'label', 'Star Size', 'label'),
            ('Sample Area', 'label'),
            ('Exptime', 'entry'),
            )

        w, b = GtkHelp.build_info(captions)
        self.w = b
        b.exptime.set_text(str(self.exptime))
        self.wdetail = b

        label = gtk.Label("Select")
        label.show()
        nb.append_page(w, label)
        nb.set_tab_reorderable(w, True)
        #nb.set_tab_detachable(w, True)

        captions = (
            ('New algorithm', 'checkbutton'),
            ('Radius', 'xlabel', '@Radius', 'spinbutton'),
            ('Threshold', 'xlabel', '@Threshold', 'entry'),
            )

        w, b = GtkHelp.build_info(captions)
        self.w.update(b)
        
        b.radius.set_tooltip_text("Radius for peak detection")
        b.threshold.set_tooltip_text("Threshold for peak detection (blank=default)")

        b.new_algorithm.set_active(self.use_new_algorithm)
        def new_alg_cb(w):
            self.use_new_algorithm = w.get_active()
        b.new_algorithm.connect('toggled', new_alg_cb)

        # radius control
        adj = b.radius.get_adjustment()
        b.radius.set_digits(2)
        b.radius.set_numeric(True)
        adj.configure(self.radius, 5.0, 200.0, 1.0, 10.0, 0)
        def chg_radius(w):
            self.radius = float(w.get_text())
            self.w.lbl_radius.set_text(str(self.radius))
            return True
        b.lbl_radius.set_text(str(self.radius))
        b.radius.connect('value-changed', chg_radius)

        # threshold control
        def chg_threshold(w):
            threshold = None
            ths = w.get_text().strip()
            if len(ths) > 0:
                threshold = float(ths)
            self.threshold = threshold
            self.w.lbl_threshold.set_text(str(self.threshold))
            return True
        b.lbl_threshold.set_text(str(self.threshold))
        b.threshold.connect('activate', chg_threshold)

        label = gtk.Label("Settings")
        label.show()
        nb.append_page(w, label)
        nb.set_tab_reorderable(w, True)
        #nb.set_tab_detachable(w, True)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(5)
        btn = gtk.Button('Ok')
        btn.connect('clicked', lambda w: self.ok())
        btns.add(btn)
        btn = gtk.Button('Cancel')
        btn.connect('clicked', lambda w: self.cancel())
        btns.add(btn)
        vbox.pack_start(btns, fill=True, expand=False)
        vbox.show_all()
        
        container.pack_start(sw, padding=0, fill=True, expand=True)

    def set_message(self, msg):
        buf = self.tw.get_buffer()
        buf.set_text(msg)
        self.tw.modify_font(self.msgFont)
            
    def instructions(self):
        self.set_message("""Please select a guide star manually.

Draw (or redraw) a region with the right mouse button.  Move the region with the left mouse button.  Press Ok or Cancel to finish.""")
            
    def start(self, future=None):
        self.callerInfo = future
        # Gather parameters
        p = future.get_data()
        
        # remove all qdas canvases
        self.withdraw_qdas_layers()
        
        # insert our canvas to fitsimage if it is not already
        try:
            obj = self.fitsimage.getObjectByTag(self.layertag)

        except KeyError:
            # Add canvas layer
            self.fitsimage.add(self.canvas, tag=self.layertag)

        self.canvas.deleteAllObjects(redraw=False)
        self.instructions()

        # Save calling parameters
        #self.w.exptime.set_text(str(p.exptime))

        self.resume()

        self.agarea = p.ag_area.upper()

        # Draw exposure region, if any
        if p.has_key('er_x1'):
            self.logger.info("Exposure region at %dx%d, %dx%d" % (
                p.er_x1, p.er_y1, p.er_x2, p.er_y2))
            self.canvas.add(CanvasTypes.CompoundObject(
                CanvasTypes.Rectangle(p.er_x1, p.er_y1, p.er_x2, p.er_y2,
                                      color=self.eregcolor),
                CanvasTypes.Text(p.er_x1, p.er_y2+4, "Exposure Range",
                                 color=self.eregcolor)),
                            redraw=False)

        tag = self.canvas.add(CanvasTypes.Rectangle(p.x1, p.y1, p.x2, p.y2,
                                                    color='cyan',
                                                    linestyle='dash'))
        self.setpickregion(self.canvas, tag)


    def pause(self):
        self.canvas.ui_setActive(False)
        
    def resume(self):
        self.canvas.ui_setActive(True)
        self.fv.showStatus("Draw a rectangle with the right mouse button")
        
    def stop(self):
        # disable canvas
        self.canvas.ui_setActive(False)

    def close(self):
        chname = self.fv.get_channelName(self.fitsimage)
        self.fv.stop_operation_channel(chname, str(self))
        return True
        
    def release_caller(self):
        try:
            self.close()
        except:
            pass
        self.callerInfo.resolve(0)
        
    def ok(self):
        self.logger.info("OK clicked.")
        p = self.callerInfo.get_data()

        try:
            obj = self.canvas.getObjectByTag(self.picktag)
        except KeyError:
            # No rectangle drawn
            # TODO: throw up a popup
            pass
            
        if obj.kind != 'compound':
            return True
        bbox  = obj.objects[0]
        point = obj.objects[1]
        p.x1, p.y1, p.x2, p.y2 = bbox.x1, bbox.y1, bbox.x2, bbox.y2

        # Save exposure time for next invocation
        try:
            self.exptime = float(self.w.exptime.get_text())
        except ValueError:
            self.exptime = 3000.0
        p.exptime = self.exptime
        
        p.result = 'ok'
        self.release_caller()

    def cancel(self):
        self.logger.info("CANCEL clicked.")
        p = self.callerInfo.get_data()
        p.result = 'cancel'

        try:
            obj = self.canvas.getObjectByTag(self.picktag)
            if obj.kind == 'compound':
                bbox  = obj.objects[0]
                point = obj.objects[1]
                label = obj.objects[2]
                bbox.color = 'red'
                label.color = 'red'
                self.canvas.redraw()
        except KeyError:
            pass
        self.release_caller()

    def show_area(self, p):
        # remove all qdas canvases
        self.withdraw_qdas_layers()
        
        # insert our canvas to fitsimage if it is not already
        try:
            obj = self.fitsimage.getObjectByTag(self.layertag)

        except KeyError:
            # Add canvas layer
            self.fitsimage.add(self.canvas, tag=self.layertag)

        self.agarea = p.ag_area.upper()
        
        self.canvas.deleteAllObjects(redraw=False)

        # Draw exposure region, if any
        if p.has_key('er_x1'):
            self.logger.info("Exposure region at %dx%d, %dx%d" % (
                p.er_x1, p.er_y1, p.er_x2, p.er_y2))
            self.canvas.add(CanvasTypes.CompoundObject(
                CanvasTypes.Rectangle(p.er_x1, p.er_y1, p.er_x2, p.er_y2,
                                      color=self.eregcolor),
                CanvasTypes.Text(p.er_x1, p.er_y2+4, "Exposure Range",
                                 color=self.eregcolor)),
                            redraw=False)

        # Draw area selection
        tag = self.canvas.add(CanvasTypes.CompoundObject(
            CanvasTypes.Rectangle(p.x1, p.y1, p.x2, p.y2,
                                  color=self.pickcolor),
            CanvasTypes.Point(p.x, p.y, 10, color='green'),
            CanvasTypes.Text(p.x1, p.y2+4, self.agarea,
                             color=self.pickcolor)),
                         redraw=True)
        self.picktag = tag

        # disable canvas
        self.stop()
        return 0
        

    def redo(self):
        obj = self.canvas.getObjectByTag(self.picktag)
        if obj.kind != 'compound':
            return True
        bbox  = obj.objects[0]
        point = obj.objects[1]
        data_x, data_y = point.x, point.y
        x1, y1, x2, y2 = bbox.x1, bbox.y1, bbox.x2, bbox.y2

        p = self.callerInfo.get_data()
        try:
            image = self.fitsimage.get_image()

            # sanity check on region
            width = bbox.x2 - bbox.x1
            height = bbox.y2 - bbox.y1
            dx = width // 2
            dy = height // 2
            if (width > self.max_len) or (height > self.max_len):
                errmsg = "Image area (%dx%d) too large!" % (
                    width, height)
                self.fv.showStatus(errmsg)
                raise Exception(errmsg)
        
            # Note: FITS coordinates are 1-based, whereas numpy FITS arrays
            # are 0-based
            fits_x, fits_y = data_x + 1, data_y + 1

            self.wdetail.sample_area.set_text('%dx%d' % (width, height))

            if self.use_new_algorithm:
                qualsize = self.iqcalc.qualsize
            else:
                qualsize = self.iqcalc.qualsize_old
            qs = qualsize(image, x1, y1, x2, y2, radius=self.radius,
                          threshold=self.threshold)
            p.x, p.y = qs.x, qs.y

            # Calculate X/Y of center of star
            obj_x = qs.objx
            obj_y = qs.objy
            fwhm = qs.fwhm
            # set region from center of detected object
            p.x1, p.y1 = max(0, p.x - dx), max(0, p.y - dy)
            width = image.width
            height = image.height
            p.x2 = min(width - 1,  p.x + dx)
            p.y2 = min(height - 1, p.y + dy)

            p.fwhm = qs.fwhm
            p.brightness = qs.brightness
            p.skylevel = qs.skylevel
            p.objx = qs.objx
            p.objy = qs.objy

            self.wdetail.fwhm.set_text('%.3f' % fwhm)
            self.wdetail.object_x.set_text('%.3f' % (obj_x+1))
            self.wdetail.object_y.set_text('%.3f' % (obj_y+1))
            self.wdetail.sky_level.set_text('%.3f' % qs.skylevel)
            self.wdetail.brightness.set_text('%.3f' % qs.brightness)

            # Mark object center on image
            point.x, point.y = int(qs.x), int(qs.y)
            point.color = 'cyan'

            # Calc RA, DEC, EQUINOX of X/Y center pixel
            ra_txt, dec_txt = image.pixtoradec(obj_x, obj_y, format='str')
            self.wdetail.ra.set_text(ra_txt)
            self.wdetail.dec.set_text(dec_txt)
            equinox = image.get_keyword('EQUINOX', 'UNKNOWN')
            self.wdetail.equinox.set_text(str(equinox))

            self.pick_qs = qs

            # TODO: Get separate FWHM for X and Y
            cdelt1, cdelt2 = image.get_keywords_list('CDELT1', 'CDELT2')
            starsize = self.iqcalc.starsize(fwhm, cdelt1, fwhm, cdelt2)
            self.wdetail.star_size.set_text('%.3f' % starsize)
            
        except Exception, e:
            point.color = 'red'
            self.logger.error("Error calculating quality metrics: %s" % (
                str(e)))
            for key in ('sky_level', 'brightness', 'star_size'):
                self.wdetail[key].set_text('')
            self.wdetail.fwhm.set_text('Failed')
            self.pick_qs = None

            # set region 
            p.x1, p.y1 = max(0, x1), max(0, y1)
            width = image.width
            height = image.height
            p.x2 = min(width - 1,  x2)
            p.y2 = min(height - 1, y2)

        self.canvas.redraw(whence=3)
        
        self.fv.showStatus("Click left mouse button to reposition pick")
        return True
    
    def update(self, canvas, action, data_x, data_y):
        try:
            obj = self.canvas.getObjectByTag(self.picktag)
            if obj.kind == 'rectangle':
                bbox = obj
            else:
                bbox  = obj.objects[0]
                point = obj.objects[1]
            self.dx = (bbox.x2 - bbox.x1) // 2
            self.dy = (bbox.y2 - bbox.y1) // 2
        except Exception, e:
            pass

        dx = self.dx
        dy = self.dy
        
        # Mark center of object and region on main image
        try:
            self.canvas.deleteObjectByTag(self.picktag, redraw=False)
        except:
            pass

        x1, y1 = data_x - dx, data_y - dy
        x2, y2 = data_x + dx, data_y + dy
        
        tag = self.canvas.add(CanvasTypes.Rectangle(x1, y1, x2, y2,
                                                    color='cyan',
                                                    linestyle='dash'),
                              redraw=False)

        self.setpickregion(self.canvas, tag)
        return True
        
    def drag(self, canvas, action, data_x, data_y):

        obj = self.canvas.getObjectByTag(self.picktag)
        if obj.kind == 'compound':
            bbox = obj.objects[0]
        elif obj.kind == 'rectangle':
            bbox = obj
        else:
            return True

        # calculate center of bbox
        wd = bbox.x2 - bbox.x1
        dw = wd // 2
        ht = bbox.y2 - bbox.y1
        dh = ht // 2
        x, y = bbox.x1 + dw, bbox.y1 + dh

        # calculate offsets of move
        dx = (data_x - x)
        dy = (data_y - y)

        # calculate new coords
        x1, y1, x2, y2 = bbox.x1+dx, bbox.y1+dy, bbox.x2+dx, bbox.y2+dy

        if (not obj) or (obj.kind == 'compound'):
            # Replace compound image with rectangle
            try:
                self.canvas.deleteObjectByTag(self.picktag, redraw=False)
            except:
                pass

            self.picktag = self.canvas.add(CanvasTypes.Rectangle(x1, y1, x2, y2,
                                                                 color='cyan',
                                                                 linestyle='dash'))
        else:
            # Update current rectangle with new coords and redraw
            bbox.x1, bbox.y1, bbox.x2, bbox.y2 = x1, y1, x2, y2
            self.canvas.redraw(whence=3)

        return True

    def setpickregion(self, canvas, tag):
        obj = canvas.getObjectByTag(tag)
        if obj.kind != 'rectangle':
            return True
        canvas.deleteObjectByTag(tag, redraw=False)

        if self.picktag:
            try:
                canvas.deleteObjectByTag(self.picktag, redraw=False)
            except:
                pass

        # determine center of rectangle
        x = obj.x1 + (obj.x2 - obj.x1) // 2
        y = obj.y1 + (obj.y2 - obj.y1) // 2
        
        tag = canvas.add(CanvasTypes.CompoundObject(
            CanvasTypes.Rectangle(obj.x1, obj.y1, obj.x2, obj.y2,
                                  color=self.pickcolor),
            CanvasTypes.Point(x, y, 10, color='red'),
            CanvasTypes.Text(obj.x1, obj.y2+4, self.agarea,
                             color=self.pickcolor)),
                         redraw=False)
        self.picktag = tag

        #self.fv.raise_tab("detail")
        self.redo()
        return True
    
    def __str__(self):
        return 'agareaselection'
    
#END
