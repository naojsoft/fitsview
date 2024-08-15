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
to and press the "Set Memo" button.  Multiple selection follows the usual
rules about holding down CTRL and/or SHIFT keys.

***Displaying an image***

Double-click on a log entry.

"""
import os
from datetime import datetime
from dateutil import tz
from collections import OrderedDict

from ginga import GingaPlugin, AstroImage
from ginga.gw import Widgets

from g2cam.INS import INSdata

__all__ = ['ObsLog']


class ObsLog(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        super(ObsLog, self).__init__(fv)

        self.chnames = []
        self.file_prefixes = []
        self.auto_scroll = True
        self.sort_hdr = 'FrameID'

        # columns to be shown in the table
        column_info = [dict(col_title="Obs Mod", fits_kwd='OBS-MOD'),
                       dict(col_title="Datatype", fits_kwd='DATA-TYP'),
                       dict(col_title=self.sort_hdr, fits_kwd='FRAMEID'),
                       dict(col_title="Object", fits_kwd='OBJECT'),
                       dict(col_title="UT", fits_kwd='UT'),
                       dict(col_title="PropId", fits_kwd='PROP-ID'),
                       dict(col_title="Exp Time", fits_kwd='EXPTIME'),
                       dict(col_title="Air Mass", fits_kwd='AIRMASS'),
                       dict(col_title="RA", fits_kwd='RA'),
                       dict(col_title="DEC", fits_kwd='DEC'),
                       dict(col_title="EQUINOX", fits_kwd='EQUINOX'),
                       dict(col_title="Memo", fits_kwd='G_MEMO'),
                       ]

        prefs = self.fv.get_preferences()
        self.settings = prefs.create_category('plugin_ObsLog')
        self.settings.add_defaults(sortable=True,
                                   color_alternate_rows=True,
                                   column_info=column_info,
                                   cache_normalized_images=True)

        self.rpt_dict = OrderedDict({})
        self.rpt_columns = []
        self.col_widths = []
        self.memo_txt = ''
        self.insconfig = INSdata()

        self.col_info = self.settings.get('column_info', [])
        # this will set rpt_columns and col_widths
        self.process_columns(self.col_info)

        self.fv.add_callback('add-image', self.incoming_data_cb)
        self.gui_up = False

    def process_columns(self, spec_lst):
        rpt_columns = []
        col_widths = []
        for dct in spec_lst:
            rpt_columns.append((dct['col_title'], dct['fits_kwd']))
            col_widths.append(dct.get('col_width', None))
        self.rpt_columns = rpt_columns
        self.col_widths = col_widths

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

        tv.setup_table(self.rpt_columns, 1, 'FRAMEID')

        # set any specified column widths
        tv.set_optimal_column_widths()
        for i, wd in enumerate(self.col_widths):
            if wd is None:
                continue
            tv.set_column_width(i, wd)

        captions = (("Memo:", 'label', "memo", 'entry', "Set Memo", 'button',
                     "Copy Memo", 'button', "Paste Memo", 'button'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w.update(b)
        vbox.add_widget(w, stretch=0)

        b.memo.set_tooltip('Set memo for selected frames')
        b.set_memo.add_callback('activated', self.set_memo_cb)
        b.set_memo.set_enabled(False)
        b.copy_memo.add_callback('activated', self.copy_memo_cb)
        #b.copy_memo.set_enabled(False)
        b.paste_memo.add_callback('activated', self.paste_memo_cb)
        #b.paste_memo.set_enabled(False)

        captions = (("Folder:", 'label', "obslog_dir", 'entry',
                     "Name:", 'label', "obslog_name", 'entryset',
                     "Type", 'combobox', "Load", 'button',
                     "merge", 'checkbutton', "Save", 'button',
                     "auto save", 'checkbutton'),
                    )
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w.update(b)
        vbox.add_widget(w, stretch=0)

        obs_log = self.settings.get('obslog_name', None)
        if obs_log is None:
            now = datetime.now(tz=tz.UTC)
            obs_log = now.strftime("obslog-%Y-%m-%d.csv")
        b.obslog_name.set_text(obs_log)
        b.obslog_name.set_tooltip('File name for observation log')
        #b.obslog_name.add_callback('activated', self.save_obslog_cb)

        b.obslog_dir.set_text("/tmp")
        b.obslog_dir.set_tooltip('Folder path for observation log')
        #b.obslog_dir.add_callback('activated', self.save_obslog_cb)

        b.type.insert_alpha("csv")
        b.type.insert_alpha("xlsx")
        b.type.set_tooltip("Format for saving/loading ObsLog")
        b.type.add_callback('activated', self.set_obslog_format_cb)

        b.load.set_tooltip("Load a saved ObsLog")
        b.load.add_callback('activated', self.load_obslog_cb)

        b.merge.set_tooltip("Merge the loaded ObsLog with current table")
        b.merge.set_state(False)

        b.save.set_tooltip("Save the ObsLog now")
        b.save.add_callback('activated', self.save_obslog_cb)

        b.auto_save.set_tooltip("Automatically save the ObsLog when new entries are added")
        b.auto_save.set_state(False)

        btns = Widgets.HBox()
        btns.set_border_width(4)
        btns.set_spacing(4)

        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        #btn.set_enabled(False)
        btns.add_widget(btn)
        btn = Widgets.Button("Save column widths")
        btn.add_callback('activated', self.record_col_widths_cb)
        btns.add_widget(btn, stretch=0)
        btn = Widgets.CheckBox("Auto scroll")
        btn.set_state(self.auto_scroll)
        btn.add_callback('activated', self.set_auto_scroll_cb)
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        vbox.add_widget(btns, stretch=0)

        container.add_widget(vbox, stretch=1)
        self.gui_up = True

        self.update_obslog()

    def replace_kwds(self, header):
        """Subclass this method to do munge the data for special reports."""
        d = dict()
        d.update(header)
        return d

    def add_to_obslog(self, header, image):
        frameid = header['FRAMEID']

        if frameid in self.rpt_dict:
            # already an entry for this frame?
            return

        # replace some kwds as needed in the table
        d = self.replace_kwds(header.asdict())

        # Hack to insure that we get the columns in the desired order
        d = OrderedDict([(kwd, d.get(kwd, ''))
                         for col, kwd in self.rpt_columns])
        self.rpt_dict[frameid] = d
        self.logger.info("adding to dict [{}]: {}".format(frameid, str(d)))

        self.update_obslog()

    def stop(self):
        self.gui_up = False

    def process_image(self, chname, header, image):
        """Override this method to do something special with the data."""
        pass

    def incoming_data_cb(self, fv, chname, image, info):
        if chname not in self.chnames:
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

        if self.auto_scroll:
            self.w.rpt_tbl.scroll_to_end()

        if self.w.auto_save.get_state():
            obslog_name = self.w.obslog_name.get_text().strip()
            if len(obslog_name) > 0:
                obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                           obslog_name)
                self.save_obslog(obslog_path)

    def save_obslog(self, filepath):
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

    def load_obslog(self, filepath, merge=False):
        try:
            import pandas as pd
        except ImportError:
            self.fv.show_error("Please install 'pandas' and "
                               "'openpyxl' to use this feature")
            return

        save_dict = dict()
        if merge:
            save_dict = dict(self.rpt_dict)

        try:
            self.logger.info("loading obslog: {}".format(filepath))

            col_hdr = [key for colname, key in self.rpt_columns]
            if filepath.endswith('.csv'):
                df = pd.read_csv(filepath, header=0, names=col_hdr,
                                 keep_default_na=False, index_col=None)

            else:
                df = pd.read_excel(filepath, header=0, names=col_hdr,
                                   index_col=None)

            res = df.to_dict('index')

            # merge loaded table (as dict) into save_dict
            for row in res.values():
                frameid = row['FRAMEID']
                d = OrderedDict([(kwd, row.get(kwd, ''))
                                 for col, kwd in self.rpt_columns])
                save_dict[frameid] = d

            rows = list(save_dict.values())
            rows.sort(key=lambda dct: dct['FRAMEID'])

            rpt_dict = OrderedDict([(row['FRAMEID'], row) for row in rows])

            self.rpt_dict = rpt_dict
            self.w.rpt_tbl.set_tree(self.rpt_dict)

        except Exception as e:
            self.logger.error("Error loading obslog: {}".format(e),
                              exc_info=True)

    def save_obslog_cb(self, w):
        obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                   self.w.obslog_name.get_text().strip())
        self.save_obslog(obslog_path)

    def load_obslog_cb(self, w):
        obslog_path = os.path.join(self.w.obslog_dir.get_text().strip(),
                                   self.w.obslog_name.get_text().strip())
        merge = self.w.merge.get_state()
        self.load_obslog(obslog_path, merge=merge)

    def get_selected(self):
        res_dict = self.w.rpt_tbl.get_selected()
        return res_dict

    def dblclick_cb(self, widget, d):
        """Switch to the image that was double-clicked in the obslog"""
        frameid = list(d.keys())[0]
        info = d[frameid]

        self.view_image(frameid, info)

    def view_image(self, frameid, info):
        chname = None
        # search for this frame in our channels; if we find it
        # switch to that image
        for _chname in self.chnames:
            if self.fv.has_channel(_chname):
                channel = self.fv.get_channel(_chname)
                if frameid in channel:
                    chname = _chname
                    channel.switch_name(frameid)
                    break

        if chname is None:
            # <-- we've apparently never loaded this file
            insname = self.insconfig.getNameByFrameId(frameid)
            filepath = os.path.join('/gen2', 'share', 'data', insname,
                                    frameid + '.fits')
            self.logger.info(f"attempting to load '{filepath}' via Gen2Int")
            pl_obj = self.fv.gpmon.get_plugin('Gen2Int')
            bnch = pl_obj.open_fits(filepath, frameid=frameid,
                                    display_image=True, wait=True)
            chname = bnch.chname

        # switch channels if we are not showing the desired one
        _channel = self.fv.get_current_channel()
        if _channel.name != chname:
            self.fv.change_channel(chname)

    def select_cb(self, widget, d):
        res = self.get_selected()
        if len(res) == 0:
            self.w.set_memo.set_enabled(False)
        else:
            self.w.set_memo.set_enabled(True)
            # grab the first memo item we can find to populate the memo box
            for frame_id, d in res.items():
                memo_txt = d['G_MEMO']
                self.w.memo.set_text(memo_txt)
                break

    def set_memo_cb(self, widget):
        memo_txt = self.w.memo.get_text().strip()
        res = self.get_selected()
        if len(res) == 0:
            self.fv.show_error("No frames selected for memo!")
            return

        for key in res.keys():
            self.rpt_dict[key]['G_MEMO'] = memo_txt

        self.update_obslog()

    def copy_memo_cb(self, widget):
        self.memo_txt = self.w.memo.get_text().strip()
        if len(self.memo_txt) > 0:
            self.fv.show_status("Memo copied")
            return

    def paste_memo_cb(self, widget):
        self.w.memo.set_text(self.memo_txt)

    def record_col_widths_cb(self, widget):
        widths = self.w.rpt_tbl.get_column_widths()
        column_info = self.settings.get('column_info')
        for i, wd in enumerate(widths):
            column_info[i]['col_width'] = wd

        # save settings
        self.settings.save()

    def set_obslog_format_cb(self, w, idx):
        ext = w.get_text()
        obslog_name = self.w.obslog_name.get_text().strip()
        name, old_ext = os.path.splitext(obslog_name)

        self.w.obslog_name.set_text(name + '.' + ext)
        #self.save_obslog_cb(None)

    def set_auto_scroll_cb(self, widget, tf):
        self.auto_scroll = tf

    def close(self):
        self.fv.stop_global_plugin(str(self))
        return True

    def __str__(self):
        return 'obslog'
