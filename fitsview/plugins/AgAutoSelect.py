#
# AgAutoSelect.py -- A VGW plugin for fits viewer
#
# E. Jeschke
#
from ginga.misc import Bunch
from ginga.rv.plugins import Catalogs

# $PYHOME imports
import astro.radec as radec

from fitsview.util import g2catalog


class AgAutoSelect(Catalogs.Catalogs):

    def __init__(self, fv, fitsimage):
        super(AgAutoSelect, self).__init__(fv, fitsimage)

        # thickness of the marked rings for FOV
        self.ring_thickness = 2
        self.colors = Bunch.Bunch(inst='magenta', outer='red', inner='red',
                                  vignette='green', probe='cyan')
        self.probe_vignette_radius = None

    def build_gui(self, container, future=None):
        super(AgAutoSelect, self).build_gui(container, future=future)

        # add blocklist feature
        self.table.add_operation("add to blocklist", self.add_blocklist)
        self.table.add_operation("rm from blocklist", self.rem_blocklist)
        self.table.btn['oprn'].append_text("add to blocklist")
        self.table.btn['oprn'].append_text("rm from blocklist")
        self.table.btn['oprn'].set_index(0)

    def start(self, future):
        self.callerInfo = future
        super(AgAutoSelect, self).start()

    def get_canvas(self):
        return self.canvas

    def plot(self, future, plotObj):
        self.callerInfo = future
        # Gather parameters
        p = future.get_data()

        self.clear_all()
        self.reset()

        # Draw the graphics for this particular foci or instrument
        plotObj.draw(self)

        # If there is a starlist waiting to be plotted, do it
        if p.starlist:
            if (self.limit_stars_to_area and
                hasattr(plotObj, 'filter_results')):
                starlist = plotObj.filter_results(p.starlist)
            else:
                starlist = p.starlist
            p.starlist = starlist
        else:
            p.starlist = []

        self.probe_vignette_radius = p.get('probe_vignette_radius', None)

        # Update GUI
        self.update_catalog(p.starlist, p.info)

        # Select top star
        if len(p.starlist) > 0:
            self.table.show_selection(p.starlist[0])

        #self.fv.update_pending(timeout=0.25)

    def highlight_object(self, obj, tag, color, redraw=True):
        x = obj.objects[0].x
        y = obj.objects[0].y
        delta = 10
        radius = obj.objects[0].radius + delta

        hilite = self.dc.CompoundObject()
        # TODO: we have to add this to the canvas first--fix this
        self.hilite.add_object(hilite)

        hilite.add_object(self.dc.Circle(x, y, radius,
                                         linewidth=4, color=color))
        # TODO: consider calling back into the plotObj for a custom
        # highlight
        if self.probe_vignette_radius is not None:
            hilite.add_object(self.dc.Circle(x, y, self.probe_vignette_radius,
                                             linewidth=2, color='green',
                                             linestyle='dash'))
        if redraw:
            self.canvas.update_canvas()


    def release_caller(self):
        self.callerInfo.resolve(0)

    def close(self):
        self.ok()

        chname = self.fv.get_channel_name(self.fitsimage)
        self.fv.stop_local_plugin(chname, str(self))
        return True

    def ok(self):
        self.logger.info("OK clicked.")
        p = self.callerInfo.get_data()

        p.result = 'ok'
        selected = self.table.get_selected()
        if len(selected) == 0:
            self.fv.show_error("No object selected.  Please select an object!")
            return False
        p.selected = selected
        self.logger.debug("returning %s" % str(p.selected))
        self.release_caller()
        return True

    def cancel(self):
        self.logger.info("CANCEL clicked.")
        p = self.callerInfo.get_data()
        p.result = 'cancel'
        self.release_caller()
        return True

    def add_blocklist(self, selected):
        self.logger.info("selected=%s" % (str(selected)))
        star = selected[0]
        g2catalog.blocklist.add_blocklist(star, self.logger)

    def rem_blocklist(self, selected):
        self.logger.info("selected=%s" % (str(selected)))
        star = selected[0]
        g2catalog.blocklist.remove_blocklist(star,  self.logger)

    def __str__(self):
        return 'agautoselect'

#END
