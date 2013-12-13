#
# Ana_UserInput.py -- ANA user input plugin for fits viewer
# 
# Eric Jeschke (eric@naoj.org)
#
import gtk
import pango
from ginga.gtkw import GtkHelp

from ginga import GingaPlugin

class Ana_UserInput(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(Ana_UserInput, self).__init__(fv, fitsimage)

    def build_gui(self, container, future=None):
        vbox1 = gtk.VBox()
        vbox = gtk.VBox()
        
        fr = gtk.Frame()
        fr.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
        fr.set_label_align(0.1, 0.5)
        fr.add(vbox)
        
        self.lbl = gtk.Label()
        vbox.pack_start(self.lbl, padding=4, fill=True, expand=True)

        vbox2 = gtk.VBox()
        self.entries = vbox2
        vbox.pack_start(vbox2, padding=4, fill=True, expand=False)

        vbox1.pack_start(fr, padding=4, fill=True, expand=True)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(3)

        btn = gtk.Button("Ok")
        btn.connect('clicked', lambda w: self.ok())
        btns.add(btn)
        btn = gtk.Button("Cancel")
        btn.connect('clicked', lambda w: self.cancel())
        btns.add(btn)
        vbox1.pack_start(btns, padding=4, fill=True, expand=False)

        vbox1.show_all()
        container.pack_start(vbox1, padding=0, fill=True, expand=False)

    def close(self):
        chname = self.fv.get_channelName(self.fitsimage)
        self.fv.stop_operation_channel(chname, str(self))
        return True
        
    def start(self, future=None):
        self.callerInfo = future
        # Gather parameters
        p = future.get_data()

        self.lbl.set_text(p.title)

        # Remove previous entries
        for w in self.entries.get_children():
            self.entries.remove(w)

        p.resDict = {}

        tbl = gtk.Table(rows=len(p.itemlist), columns=2)
        tbl.set_row_spacings(2)
        tbl.set_col_spacings(2)

        row = 0
        for name, val in p.itemlist:
            lbl = gtk.Label(name)
            lbl.set_alignment(1.0, 0.5)
            ent = gtk.Entry()
            val_s = str(val)
            ent.set_text(val_s)
            p.resDict[name] = ent

            tbl.attach(lbl, 0, 1, row, row+1, xoptions=gtk.FILL)
            tbl.attach(ent, 1, 2, row, row+1, xoptions=gtk.EXPAND|gtk.FILL)
            row += 1

        tbl.show_all()
        self.entries.pack_start(tbl, fill=True, expand=False, padding=2)

    def resume(self):
        pass
    
    def release_caller(self):
        try:
            self.close()
        except:
            pass
        self.callerInfo.resolve(0)
        
    def ok(self):
        p = self.callerInfo.get_data()

        p.result = 'ok'
        # Read out the entry widgets
        d = {}
        for key, ent in p.resDict.items():
            s = ent.get_text()
            d[key] = s

        p.resDict = d
        self.logger.info("OK clicked, items=%s" % (d))
        self.release_caller()

    def cancel(self):
        p = self.callerInfo.get_data()
        p.result = 'cancel'
        p.resDict = {}
        self.logger.info("CANCEL clicked.")
        self.release_caller()

    def stop(self):
        pass
    
    def redo(self):
        pass
    
    def __str__(self):
        return 'ana_userinput'
    
#END
