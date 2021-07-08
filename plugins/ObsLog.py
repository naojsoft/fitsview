# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
"""Observation Log plugin.

**Plugin Type: Global**

``ObsLog`` is a global plugin. Only one instance can be opened.

**Usage**

***Saving the log to a file***

Put in values for the Observation Log folder and filename.  The format
of the file saved will depend on the file extension of the filename;
use the type selector combobox to pick the right extension:

* csv:
* xlsx: MS Excel file format

The file is rewritten out every time a new entry is added to the log

***Adding a memo to one or more log entries***

Write a memo in the memo box.  Select one or more frames to add the memo
to and press the "Add Memo" button.  Multiple selection follows the usual
rules about holding down CTRL and/or SHIFT keys.

***Displaying an image***

Double-click on a log entry.

"""
import os
from collections import OrderedDict

from ginga import GingaPlugin, AstroImage
from ginga.gw import Widgets

__all__ = ['ObsLog']


class ObsLog(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        super(ObsLog, self).__init__(fv)

        self.chname = None
        self.file_prefixes = []

        # columns to be shown in the table
        columns = [("Obs Mod", 'OBS-MOD'),
                   ("Datatype", 'DATA-TYP'),
                   ("FrameID", 'FRAMEID'),
                   ("Object", 'OBJECT'),
                   ("UT", 'UT'),
                   ("PropId", 'PROP-ID'),
                   ("Exp Time", 'EXPTIME'),
                   ("Air Mass", 'AIRMASS'),
                   #("Pos Ang", 'INST-PA'),
                   #("Ins Rot", 'INSROT'),
                   #("Foc Val", 'FOC-VAL'),
                   #("Filter01", 'FILTER01'),
                   #("Filter02", 'FILTER02'),
                   #("Filter03", 'FILTER03'),
                   ("RA", 'RA'),
                   ("DEC", 'DEC'),
                   ("EQUINOX", 'EQUINOX'),
                   ("Memo", 'G_MEMO'),
                  ]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_ObsLog')
        self.settings.add_defaults(sortable=True,
                                   color_alternate_rows=True,
                                   #max_rows_for_col_resize=5000,
                                   report_columns=columns,
                                   cache_normalized_images=True)
        self.settings.load(onError='silent')

        self.rpt_dict = OrderedDict({})
        self.rpt_columns = []

        self.fv.add_callback('add-image', self.incoming_data_cb)
        self.gui_up = False

    def build_gui(self, container):
        vbox = Widgets.VBox()
        vbox.set_border_width(1)
        vbox.set_spacing(1)

        tv = Widgets.TreeView(sortable=self.settings.get('sortable'),
                              use_alt_row_color=self.settings.get('color_alternate_rows'),
                              selection='multiple')
        self.w.rpt_tbl = tv
        vbox.add_widget(tv, stretch=1)

        tv.add_callback('activated', self.dblclick_cb)
        tv.add_callback('selected', self.select_cb)

        self.rpt_columns = self.settings.get('report_columns')
        tv.setup_table(self.rpt_columns, 1, 'FRAMEID')

        captions = (("Memo:", 'label', "memo", 'entry', 'Add Memo', 'button'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w.update(b)
        vbox.add_widget(w, stretch=0)

        b.memo.set_tooltip('Set memo for selected frames')
        b.add_memo.add_callback('activated', self.add_memo_cb)
        b.add_memo.set_enabled(False)

        captions = (("Folder:", 'label', "obslog_dir", 'entry',
                     "Name:", 'label', "obslog_name", 'entryset',
                     "Type", 'combobox', "Load", 'button'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w.update(b)
        vbox.add_widget(w, stretch=0)

        obs_log = self.settings.get('obslog_name', None)
        if obs_log is None:
            obs_log = ''
        b.obslog_name.set_text(obs_log)
        b.obslog_name.set_tooltip('File name for observation log')
        b.obslog_name.add_callback('activated', self.write_obslog_cb)

        b.obslog_dir.set_text("/tmp")
        b.obslog_dir.set_tooltip('Folder path for observation log')
        b.obslog_dir.add_callback('activated', self.write_obslog_cb)

        b.type.insert_alpha("csv")
        b.type.insert_alpha("xlsx")
        b.type.set_tooltip("Format for saving/loading ObsLog")
        b.type.add_callback('activated', self.set_obslog_format_cb)

        b.load.set_tooltip("Load a saved ObsLog")
        b.load.add_callback('activated', self.load_obslog_cb)

        btns = Widgets.HBox()
        btns.set_border_width(4)
        btns.set_spacing(4)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btn.set_enabled(False)
        btns.add_widget(btn)
        btn = Widgets.Button("Help")
        btn.add_callback('activated', lambda w: self.help())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        vbox.add_widget(btns, stretch=0)

        container.add_widget(vbox, stretch=1)
        self.gui_up = True

    def replace_kwds(self, header):
        """Subclass this method to do munge the data for special reports."""
        d = dict()
        d.update(header)
        return d

    def add_to_obslog(self, header, image):
        frameid = header['FRAMEID']

        # replace some kwds as needed in the table
        d = self.replace_kwds(header.asdict())

        # Hack to insure that we get the columns in the desired order
        d = OrderedDict([(kwd, d.get(kwd, ''))
                         for col, kwd in self.rpt_columns])
        self.rpt_dict[frameid] = d

        self.update_obslog()

    def stop(self):
        self.gui_up = False

    def process_image(self, chname, header, image):
        """Override this method to do something special with the data."""
        pass

    def incoming_data_cb(self, fv, chname, image, info):
        if chname != self.chname:
            return

        imname = image.get('name', None)
        if imname is None:
            return

        # only accepted list of frames
        accepted = False
        for prefix in self.file_prefixes:
            if imname.startswith(prefix):
                accepted = True
                break
        if not accepted:
            return

        header = image.get_header()

        # add image to obslog
        self.fv.gui_do(self.add_to_obslog, header, image)

        try:
            self.process_image(chname, header, image)

        except Exception as e:
            self.logger.error("Failed to process image: {}".format(e),
                              exc_info=True)

    def update_obslog(self):
        if not self.gui_up:
            return

        self.w.rpt_tbl.set_tree(self.rpt_dict)

        obslog_name = self.w.obslog_name.get_text().strip()
        if len(obslog_name) > 0:
            obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                       obslog_name)
            self.write_obslog(obslog_path)

    def write_obslog(self, filepath):
        if len(self.rpt_dict) == 0:
            return

        try:
            import pandas as pd
        except ImportError:
            self.fv.show_error("Please install 'pandas' and "
                               "'openpyxl' to use this feature")
            return

        try:
            self.logger.info("writing obslog: {}".format(filepath))

            col_hdr = [colname for colname, key in self.rpt_columns]
            rows = [list(d.values()) for d in self.rpt_dict.values()]
            df = pd.DataFrame(rows, columns=col_hdr)

            if filepath.endswith('.csv'):
                df.to_csv(filepath, index=False, header=True)

            else:
                df.to_excel(filepath, index=False, header=True)

        except Exception as e:
            self.logger.error("Error writing obslog: {}".format(e),
                              exc_info=True)

    def load_obslog(self, filepath):
        try:
            import pandas as pd
        except ImportError:
            self.fv.show_error("Please install 'pandas' and "
                               "'openpyxl' to use this feature")
            return

        try:
            self.logger.info("loading obslog: {}".format(filepath))

            col_hdr = [colname for colname, key in self.rpt_columns]
            if filepath.endswith('.csv'):
                df = pd.read_csv(filepath, header=0, #names=col_hdr,
                                   index_col=None)

            else:
                df = pd.read_excel(filepath, header=0, #names=col_hdr,
                                   index_col=None)

            self.rpt_dict = OrderedDict({})

            res = df.to_dict('index')

            for row in res.values():
                frameid = row['FrameID']
                d = OrderedDict([(kwd, row.get(col, ''))
                                 for col, kwd in self.rpt_columns])
                self.rpt_dict[frameid] = d

            self.w.rpt_tbl.set_tree(self.rpt_dict)

        except Exception as e:
            self.logger.error("Error loading obslog: {}".format(e),
                              exc_info=True)

    def write_obslog_cb(self, w):
        obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                   self.w.obslog_name.get_text().strip())
        self.write_obslog(obslog_path)

    def load_obslog_cb(self, w):
        obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                   self.w.obslog_name.get_text().strip())
        self.load_obslog(obslog_path)

    def get_selected(self):
        res_dict = self.w.rpt_tbl.get_selected()
        return res_dict

    def dblclick_cb(self, widget, d):
        """Switch to the image that was double-clicked in the obslog"""
        frameid = list(d.keys())[0]
        info = d[frameid]

        self.view_image(frameid, info)

    def view_image(self, imname, info):
        chname = self.chname
        channel = self.fv.get_current_channel()
        if channel.name != chname:
            channel = self.fv.get_channel(chname)
            self.fv.change_channel(chname)

        channel.switch_name(imname)

    def select_cb(self, widget, d):
        res = self.get_selected()
        if len(res) == 0:
            self.w.add_memo.set_enabled(False)
        else:
            self.w.add_memo.set_enabled(True)

    def add_memo_cb(self, widget):
        memo_txt = self.w.memo.get_text().strip()
        res = self.get_selected()
        if len(res) == 0:
            self.fv.show_error("No frames selected for memo!")
            return

        for key in res.keys():
            self.rpt_dict[key]['G_MEMO'] = memo_txt

        self.w.rpt_tbl.set_tree(self.rpt_dict)

    def set_obslog_format_cb(self, w, idx):
        ext = w.get_text()
        obslog_name = self.w.obslog_name.get_text().strip()
        name, old_ext = os.path.splitext(obslog_name)

        self.w.obslog_name.set_text(name + '.' + ext)
        self.write_obslog_cb(None)

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def __str__(self):
        return 'obslog'
