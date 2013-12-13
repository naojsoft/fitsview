#
# Ana_Confirmation.py -- ANA Confirmation plugin for fits viewer
# 
# Eric Jeschke (eric@naoj.org)
#
import gtk
from ginga.gtkw import GtkHelp

from ginga import GingaPlugin

class Ana_Confirmation(GingaPlugin.LocalPlugin):

    def __init__(self, fv, fitsimage):
        # superclass defines some variables for us, like logger
        super(Ana_Confirmation, self).__init__(fv, fitsimage)

    def build_gui(self, container, future=None):
        vbox1 = gtk.VBox()
        vbox = gtk.VBox()
        
        fr = gtk.Frame()
        fr.set_shadow_type(gtk.SHADOW_ETCHED_OUT)
        fr.set_label_align(0.1, 0.5)
        fr.add(vbox)
        
        self.lbl = gtk.Label()
        vbox.pack_start(self.lbl, padding=4, fill=True, expand=True)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(3)
        self.btns = btns
        vbox.pack_start(btns, padding=4, fill=True, expand=False)

        vbox1.pack_start(fr, padding=4, fill=True, expand=True)

        btns = gtk.HButtonBox()
        btns.set_layout(gtk.BUTTONBOX_START)
        btns.set_spacing(3)

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

        # Remove previous buttons
        for w in self.btns.get_children():
            self.btns.remove(w)

        def cb(index):
            return lambda w: self.ok(index)

        items = p.dialog.split()
        index = 1
        for name in items:
            btn = gtk.Button(name)
            btn.connect('clicked', cb(index))
            self.btns.add(btn)
            index += 1

    def resume(self):
        pass
    
    def release_caller(self):
        try:
            self.close()
        except:
            pass
        self.callerInfo.resolve(0)
        
    def ok(self, index):
        self.logger.info("OK clicked, index=%d" % (index))
        p = self.callerInfo.get_data()

        p.result = 'ok'
        p.index  = index
        self.release_caller()

    def cancel(self):
        self.logger.info("CANCEL clicked.")
        p = self.callerInfo.get_data()
        p.result = 'cancel'
        self.release_caller()

    def stop(self):
        pass
    
    def redo(self):
        pass
    
    def __str__(self):
        return 'ana_confirmation'
    
#END
