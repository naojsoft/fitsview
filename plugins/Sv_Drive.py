#
# Sv_Drive.py -- Object/destination calculation plugin for fits viewer
# 
#[ Eric Jeschke (eric@naoj.org) --
#  Last edit: Tue Sep  9 12:01:04 HST 2014
#]
#
import gtk
from ginga.misc import Bunch

from ginga.gtkw import GtkHelp
from ginga.gtkw import ImageViewCanvasGtk
from ginga.gtkw import ImageViewCanvasTypesGtk as CanvasTypes
from ginga import GingaPlugin

# Local application imports
from util import g2calc


class QDASPlugin(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        super(QDASPlugin, self).__init__(fv, fitsimage)

        self.layertag = 'qdas-svdrive'

        canvas = CanvasTypes.DrawingCanvas()
        canvas.enable_draw(True)
        self.canvas = canvas

        #self.canvas.set_callback('button-press', self.update)
        canvas.set_callback('key-press', self.keydown)
        canvas.set_callback('cursor-down', self.btndown)
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
        


class Sv_Drive(QDASPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(Sv_Drive, self).__init__(fv, fitsimage)

        self.objtag = None
        self.dsttag = None
        self.regiontag = None

        self.dst_x = 0
        self.dst_y = 0
        self.obj_x = 0
        self.obj_y = 0
        self.x1 = 0
        self.y1 = 0
        self.x2 = 0
        self.y2 = 0
        self.recenter = False
        self.width = None
        self.height = None
        self.isDst = True

        self.use_new_algorithm = False
        self.radius = 10
        self.threshold = None
        self.min_fwhm = 2.0
        self.max_fwhm = 50.0
        self.min_ellipse = 0.5
        self.edgew = 0.01

        # this is the maximum size a side can be in bounding box
        self.max_len = 1024

        # For image FWHM type calculations
        self.iqcalc = g2calc.IQCalc(self.logger)

        self.have_gui = False

    def get_dst(self):
        return (self.dst_x, self.dst_y)
    
    def get_obj(self):
        return (self.obj_x, self.obj_y)
    
    def get_region(self):
        return (self.x1, self.y1, self.x2, self.y2)
    
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
            ('Dst X', 'entry', 'Dst Y', 'entry'),
            ('Obj X', 'entry', 'Obj Y', 'entry'),
            ('X1', 'entry', 'Y1', 'entry'),
            ('X2', 'entry', 'Y2', 'entry'),
            ('Update Pos', 'button', 'Recenter', 'checkbutton'),
            ('Frame', 'combobox'),
            ('Dst RA', 'label', 'Dst DEC', 'label'),
            ('Obj RA', 'label', 'Obj DEC', 'label'),
            #('Sky Level', 'label', 'Brightness', 'label'), 
            #('FWHM', 'label', 'Star Size', 'label'),
            )

        w, b = GtkHelp.build_info(captions)
        self.w = b
        self.w.update_pos.connect('clicked', lambda w: self.update_positions_cb())
        self.w.dst_x.connect('activate', lambda w: self.update_positions_cb())
        self.w.dst_y.connect('activate', lambda w: self.update_positions_cb())
        self.w.obj_x.connect('activate', lambda w: self.update_positions_cb())
        self.w.obj_y.connect('activate', lambda w: self.update_positions_cb())
        self.w.recenter.sconnect("toggled", self.toggle_recenter)
        self.w.frame.sconnect("changed", self.change_frame_cb)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(5)
        btn = GtkHelp.RadioButton(None, "Place Destination")
        self.w.r_dst = btn
        btns.add(btn)
        btn = GtkHelp.RadioButton(btn, "Place Object")
        self.w.r_obj = btn
        btns.add(btn)
        if self.isDst:
            self.w.r_dst.set_active(True)
        else:
            self.w.r_obj.set_active(True)
        self.w.r_dst.connect("toggled", lambda w: self.toggle_dstsrc_cb())
        w.pack_start(btns, fill=True, expand=False)
        
        label = gtk.Label("Select")
        label.show()
        nb.append_page(w, label)
        nb.set_tab_reorderable(w, True)
        #nb.set_tab_detachable(w, True)
        
        captions = (
            ('New algorithm', 'checkbutton'),
            ('Radius', 'xlabel', '@Radius', 'spinbutton'),
            ('Threshold', 'xlabel', '@Threshold', 'entry'),
            ('Min FWHM', 'xlabel', '@Min FWHM', 'spinbutton'),
            ('Max FWHM', 'xlabel', '@Max FWHM', 'spinbutton'),
            ('Ellipticity', 'xlabel', '@Ellipticity', 'entry'),
            ('Edge', 'xlabel', '@Edge', 'entry'),
            )

        w, b = GtkHelp.build_info(captions)
        self.w.update(b)

        b.radius.set_tooltip_text("Radius for peak detection")
        b.threshold.set_tooltip_text("Threshold for peak detection (blank=default)")
        b.min_fwhm.set_tooltip_text("Minimum FWHM for selection")
        b.max_fwhm.set_tooltip_text("Maximum FWHM for selection")
        b.ellipticity.set_tooltip_text("Minimum ellipticity for selection")
        b.edge.set_tooltip_text("Minimum edge distance for selection")

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

        # min fwhm
        adj = b.min_fwhm.get_adjustment()
        b.min_fwhm.set_digits(2)
        b.min_fwhm.set_numeric(True)
        adj.configure(self.min_fwhm, 0.1, 200.0, 0.1, 1, 0)
        def chg_min(w):
            self.min_fwhm = w.get_value()
            self.w.lbl_min_fwhm.set_text(str(self.min_fwhm))
            return True
        b.lbl_min_fwhm.set_text(str(self.min_fwhm))
        b.min_fwhm.connect('value-changed', chg_min)

        # max fwhm
        adj = b.max_fwhm.get_adjustment()
        b.max_fwhm.set_digits(2)
        b.max_fwhm.set_numeric(True)
        adj.configure(self.max_fwhm, 0.1, 200.0, 0.1, 1, 0)
        def chg_max(w):
            self.max_fwhm = w.get_value()
            self.w.lbl_max_fwhm.set_text(str(self.max_fwhm))
            return True
        b.lbl_max_fwhm.set_text(str(self.max_fwhm))
        b.max_fwhm.connect('value-changed', chg_max)

        # Ellipticity control
        def chg_ellipticity(w):
            minellipse = None
            val = w.get_text().strip()
            if len(val) > 0:
                minellipse = float(val)
            self.min_ellipse = minellipse
            self.w.lbl_ellipticity.set_text(str(self.min_ellipse))
            return True
        b.lbl_ellipticity.set_text(str(self.min_ellipse))
        b.ellipticity.connect('activate', chg_ellipticity)

        # Edge control
        def chg_edgew(w):
            edgew = None
            val = w.get_text().strip()
            if len(val) > 0:
                edgew = float(val)
            self.edgew = edgew
            self.w.lbl_edge.set_text(str(self.edgew))
            return True
        b.lbl_edge.set_text(str(self.edgew))
        b.edge.connect('activate', chg_edgew)

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
        
        cw = container.get_widget()
        cw.pack_start(sw, padding=0, fill=True, expand=True)
        self.have_gui = True

    def set_message(self, msg):
        buf = self.tw.get_buffer()
        buf.set_text(msg)
        self.tw.modify_font(self.msgFont)
            
    def instructions(self):
        #self.set_message("""Place Destination by clicking left mouse button.  Draw a region with the right mouse button around the Object.  Press Ok or Cancel to finish.""")
        self.set_message("""Please mark object and destination.""")
            
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
        if p.has_key('msg'):
            self.set_message(p.msg)
        else:
            self.instructions()

        self.recenter = p.get('recenter', False)
        self.w.recenter.set_active(self.recenter)
        self.width = p.get('width', None)
        self.height = p.get('height', None)

        # Change the framelist
        self.frames = p.get('framelist', [])
        model = self.w.frame.get_model()
        model.clear()
        for frameid in self.frames:
            self.w.frame.append_text(frameid)
        self.w.frame.set_active(0)
        
        try:
            # IMPORTANT: Assume all coords have been adjusted from FITS
            # or CCD coords to data coords (-1)
            if p.dst_x != None:
                self.place_dst(self.canvas, p.dst_x, p.dst_y)
            else:
                self.place_dst(self.canvas, self.dst_x, self.dst_y)

            if p.obj_x != None:
                self.place_obj(self.canvas, p.obj_x, p.obj_y)
            else:
                self.place_obj(self.canvas, self.obj_x, self.obj_y)

            if p.x1 != None:
                error = False
                if p.has_key('autoerr'):
                    error = p.autoerr
                self.place_region(self.canvas, p.x1, p.y1, p.x2, p.y2,
                                  error=error)
                
        except Exception, e:
            self.logger.error("Error placing dst and objs: %s" % (
                str(e)))
            # carry on...
            
        self.resume()

    def pause(self):
        self.canvas.ui_setActive(False)
        
    def resume(self):
        # turn off any mode user may be in
        self.modes_off()

        self.canvas.ui_setActive(True)

        
    def stop(self):
        # remove the canvas from the image
        self.canvas.ui_setActive(False)

    def close(self):
        chname = self.fv.get_channelName(self.fitsimage)
        self.fv.stop_local_plugin(chname, str(self))
        self.have_gui = False
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
            self.dst_x = float(self.w.dst_x.get_text()) - 1
            self.dst_y = float(self.w.dst_y.get_text()) - 1
            self.obj_x = float(self.w.obj_x.get_text()) - 1
            self.obj_y = float(self.w.obj_y.get_text()) - 1
            self.x1 = float(self.w.x1.get_text()) - 1
            self.y1 = float(self.w.y1.get_text()) - 1
            self.x2 = float(self.w.x2.get_text()) - 1
            self.y2 = float(self.w.y2.get_text()) - 1
            ## pt = self.canvas.getObjectByTag(self.dsttag)
            ## self.dst_x, self.dst_y = pt.objects[0].x, pt.objects[0].y
            
            ## pt = self.canvas.getObjectByTag(self.objtag)
            ## self.obj_x, self.obj_y = pt.objects[0].x, pt.objects[0].y
            
            ## rect = self.canvas.getObjectByTag(self.regiontag)
            ## (self.x1, self.y1, self.x2, self.y2) = (
            ##     rect.objects[0].x1, rect.objects[0].y1,
            ##     rect.objects[0].x2, rect.objects[0].y2)

            p.dst_x = self.dst_x
            p.dst_y = self.dst_y
            p.obj_x = self.obj_x
            p.obj_y = self.obj_y
            p.x1 = self.x1
            p.y1 = self.y1
            p.x2 = self.x2
            p.y2 = self.y2
            idx = self.w.frame.get_active()
            if idx >= 0:
                p.frameid = self.frames[idx]
            else:
                p.frameid = None
            p.result = 'ok'

        except Exception, e:
            p.result = 'error'
            p.errmsg = "Error collecting dst and obj coords: %s" % (
                str(e))
            self.logger.error(p.errmsg)
            
        self.release_caller()

    def cancel(self):
        self.logger.info("CANCEL clicked.")
        p = self.callerInfo.get_data()
        p.result = 'cancel'

        self.release_caller()

    def toggle_recenter(self, w):
        recenter = w.get_active()
        self.recenter = recenter

    def redo(self):
        pass
    
    
    def place_dst(self, canvas, data_x, data_y):
        if self.dsttag:
            try:
                canvas.deleteObjectByTag(self.dsttag, redraw=False)
            except:
                pass

        x, y = data_x, data_y
        
        self.dsttag = canvas.add(CanvasTypes.CompoundObject(
            CanvasTypes.Point(x, y, 10, color='green'),
            CanvasTypes.Text(x+4, y, "Dst",
                             color='green')),
                                 redraw=False)
        
        canvas.redraw(whence=3)
        if self.have_gui:
            self.record_dst(data_x, data_y)
        
    def record_dst(self, data_x, data_y):
        self.w.dst_x.set_text('%.3f' % (data_x+1))
        self.w.dst_y.set_text('%.3f' % (data_y+1))

        # Calc RA, DEC, EQUINOX of X/Y dst pixel
        image = self.fitsimage.get_image()
        try:
            ra_txt, dec_txt = image.pixtoradec(data_x, data_y, format='str')
        except Exception, e:
            self.logger.error("Error calculating ra/dec of dst: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'
            
        self.w.dst_ra.set_text(ra_txt)
        self.w.dst_dec.set_text(dec_txt)
            
    def place_obj(self, canvas, data_x, data_y):
        if self.objtag:
            try:
                canvas.deleteObjectByTag(self.objtag, redraw=False)
            except:
                pass

        x, y = data_x, data_y
        
        # Mark object center on image
        self.objtag = canvas.add(CanvasTypes.CompoundObject(
            CanvasTypes.Point(x, y, 10, color='cyan'),
            CanvasTypes.Text(x+4, y, "Object",
                             color='green')),
                                 redraw=False)

        canvas.redraw(whence=3)
        if self.have_gui:
            self.record_obj(data_x, data_y)
        
    def record_obj(self, data_x, data_y):
        self.w.obj_x.set_text('%.3f' % (data_x+1))
        self.w.obj_y.set_text('%.3f' % (data_y+1))

        # Calc RA, DEC, EQUINOX of X/Y dst pixel
        image = self.fitsimage.get_image()
        try:
            ra_txt, dec_txt = image.pixtoradec(data_x, data_y, format='str')
        except Exception, e:
            self.logger.error("Error calculating ra/dec of object: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'

        self.w.obj_ra.set_text(ra_txt)
        self.w.obj_dec.set_text(dec_txt)
            

    def place_region(self, canvas, x1, y1, x2, y2, error=False):
        if self.regiontag:
            try:
                canvas.deleteObjectByTag(self.regiontag, redraw=False)
            except:
                pass

        color = 'cyan'; style = 'solid'
        if error:
            color = 'red'; style = 'dash'
        # Mark acquisition region on image
        self.regiontag = canvas.add(CanvasTypes.CompoundObject(
            CanvasTypes.Rectangle(x1, y1, x2, y2, color=color,
                                  linestyle=style),
            CanvasTypes.Text(x1, y2+4, "Target Acquisition",
                             color=color)),
                                 redraw=False)

        canvas.redraw(whence=3)
        if self.have_gui:
            self.record_region(x1, y1, x2, y2)

    def record_region(self, x1, y1, x2, y2):
        self.w.x1.set_text('%.3f' % (x1+1))
        self.w.y1.set_text('%.3f' % (y1+1))
        self.w.x2.set_text('%.3f' % (x2+1))
        self.w.y2.set_text('%.3f' % (y2+1))

    def update_positions_cb(self):
        dst_x = float(self.w.dst_x.get_text()) - 1
        dst_y = float(self.w.dst_y.get_text()) - 1
        obj_x = float(self.w.obj_x.get_text()) - 1
        obj_y = float(self.w.obj_y.get_text()) - 1
        x1 = float(self.w.x1.get_text()) - 1
        y1 = float(self.w.y1.get_text()) - 1
        x2 = float(self.w.x2.get_text()) - 1
        y2 = float(self.w.y2.get_text()) - 1

        self.place_dst(self.canvas, dst_x, dst_y)
        self.place_obj(self.canvas, obj_x, obj_y)
        self.place_region(self.canvas, x1, y1, x2, y2)
        return True

    def toggle_dstsrc_cb(self):
        self.isDst = self.w.r_dst.get_active()
        
    def btndown(self, canvas, action, data_x, data_y):
        self.logger.debug("Setting mark at %d,%d isDst=%s" % (
            data_x, data_y, self.isDst))
        if self.isDst:
            self.place_dst(canvas, data_x, data_y)
        else:
            self.place_obj(canvas, data_x, data_y)
        return True
        

    def keydown(self, canvas, keyname):
        if keyname == 'space':
            # toggle dst and obj placement
            isDst = self.w.r_dst.get_active()
            if isDst:
                self.w.r_obj.set_active(True)
            else:
                self.w.r_dst.set_active(True)
            return True
        elif keyname == 's':
            data_x, data_y = self.fitsimage.get_last_data_xy()
            # TODO: do something more sophisticated, like sample a region
            threshold = self.fitsimage.get_data(data_x, data_y)
            self.threshold = threshold
            self.w.lbl_threshold.set_text(str(self.threshold))
            return True

    def setpickregion(self, canvas, tag):
        bbox = canvas.getObjectByTag(tag)
        if bbox.kind != 'rectangle':
            return True
        canvas.deleteObjectByTag(tag, redraw=False)

        result = self.check_region(bbox.x1, bbox.y1, bbox.x2, bbox.y2,
                                   width=self.width, height=self.height,
                                   recenter=self.recenter)
        return True

    def check_region(self, x1, y1, x2, y2, width=None, height=None,
                     recenter=False):
        canvas = self.canvas
        if self.objtag:
            try:
                canvas.deleteObjectByTag(self.objtag, redraw=False)
            except:
                pass

        # sanity check on region
        if not width:
            width = x2 - x1
        if not height:
            height = y2 - y1
        dx = width // 2
        dy = height // 2

        result = False
        try:
            if (width > self.max_len) or (height > self.max_len):
                errmsg = "Image area (%dx%d) too large!" % (
                    width, height)
                self.set_message(errmsg)
                raise Exception(errmsg)

            image = self.fitsimage.get_image()

            if self.use_new_algorithm:
                qs = self.iqcalc.qualsize(image, x1, y1, x2, y2,
                                          radius=self.radius,
                                          threshold=self.threshold,
                                          minfwhm=self.min_fwhm,
                                          maxfwhm=self.max_fwhm,
                                          minelipse=self.min_ellipse,
                                          edgew=self.edgew)
            else:
                qs = self.iqcalc.qualsize_old(image, x1, y1, x2, y2)

            # Calculate X/Y of center of star
            obj_x = qs.objx
            obj_y = qs.objy
            self.logger.info("object center is x,y=%f,%f" % (obj_x, obj_y))
            #fwhm = qs.fwhm

            # Mark object center on image
            self.objtag = canvas.add(CanvasTypes.CompoundObject(
                CanvasTypes.Point(obj_x, obj_y, 10, color='cyan'),
                CanvasTypes.Text(obj_x+4, obj_y, "Object",
                                 color='green')),
                ## CanvasTypes.Rectangle(x1, y1, x2, y2, color='cyan')),
                                     redraw=False)

            if recenter:
                x1, y1 = max(0, obj_x - dx), max(0, obj_y - dy)
                x2 = min(image.width-1,  obj_x + dx)
                y2 = min(image.height-1, obj_y + dy)
                
            self.place_region(self.canvas, x1, y1, x2, y2)

            self.set_message("Automatic target reacquisition succeeded.")
            result = True

        except Exception, e:
            self.logger.error("Error calculating quality metrics: %s" % (
                str(e)))
            # determine center of rectangle
            obj_x = x1 + dx
            obj_y = y1 + dy

            self.objtag = canvas.add(CanvasTypes.CompoundObject(
                CanvasTypes.Point(obj_x, obj_y, 10, color='red'),
                CanvasTypes.Text(obj_x+4, obj_y, "Object",
                                 color='green')),
                ## CanvasTypes.Rectangle(x1, y1, x2, y2, color='cyan',
                ##                       linestyle='dash')),
                                     redraw=False)
            self.place_region(self.canvas, x1, y1, x2, y2, error=True)
            self.set_message("Automatic target reacquisition failed: %s" % (
                str(e)))

        self.w.obj_x.set_text('%.3f' % (obj_x+1))
        self.w.obj_y.set_text('%.3f' % (obj_y+1))

        # Calc RA, DEC, EQUINOX of X/Y object pixel
        try:
            ra_txt, dec_txt = image.pixtoradec(obj_x, obj_y, format='str')
        except Exception, e:
            self.logger.error("Error calculating ra/dec of obj: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'
            result = False
        self.w.obj_ra.set_text(ra_txt)
        self.w.obj_dec.set_text(dec_txt)
            
        canvas.redraw(whence=3)
        # Set pan position to selected object
        #self.fitsimage.panset_xy(obj_x, obj_y, redraw=False)
        return result


    def acquire_region(self, p):
        """This method is called if automatic region selection succeeded.
        """
        # remove all qdas canvases
        self.withdraw_qdas_layers()
        
        # insert our canvas to fitsimage if it is not already
        try:
            obj = self.fitsimage.getObjectByTag(self.layertag)

        except KeyError:
            # Add canvas layer
            self.fitsimage.add(self.canvas, tag=self.layertag)

        self.canvas.deleteAllObjects(redraw=False)

        print "placing dst"
        self.place_dst(self.canvas, p.dst_x, p.dst_y)

        print "placing obj"
        self.place_obj(self.canvas, p.obj_x, p.obj_y)

        print "placing region"
        self.place_region(self.canvas, p.x1, p.y1, p.x2, p.y2)
        
        # Set pan position to object
        #self.fitsimage.panset_xy(p.obj_x, p.obj_y, redraw=False)

        # disable canvas if in automatic mode
        if p.mode in ('auto', 'semiauto'):
            #self.set_message("Automatic target reacquisition succeeded.")
            self.fv.showStatus("Automatic target reacquisition succeeded.")
            #self.stop()
            self.close()
        elif p.mode == 'override':
            p.msg = "Automatic target reacquisition succeeded. Please confirm target reacquisition."
        return 0
        
    def change_frame_cb(self, w):
        p = self.callerInfo.get_data()
        index = w.get_active()
        frameid = self.frames[index]
        self.logger.debug("changing frame to '%s'" % frameid)
        p.image = p.load_frame(p.src_channel, frameid, p.dst_channel)

    def __str__(self):
        return 'sv_drive'
    
#END
