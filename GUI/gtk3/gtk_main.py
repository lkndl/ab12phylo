import sys

import gi, re
from pathlib import Path
from argparse import Namespace
import logging

# from GUI.gtk3.dataset import ab12phylo_dataset

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk, Gdk as gdk

# BASE_DIR = Path('/home/quirin/PYTHON/AB12PHYLO/GUI/gtk3/gtk_main.py').resolve().parents[2]
BASE_DIR = Path(__file__).resolve().parents[2]
__verbose__, __info__ = 1, 0

# set the icon theme
gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


class gui(gtk.Window):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'
    FILETYPES = ['.ab1', '.seq', '.fasta']

    def __init__(self, data):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(900, 600)
        self.set_size_request(640, 480)
        self.log.debug('GTK Window initialized')

        self.data = data

        # fetch all named objects from the .glade XML
        widgets = dict()
        for widget in gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                widgets[widget.get_name()] = widget
        self.interface = Namespace(**widgets)
        self.log.debug('Fetched control elements')

        # fetch the notebook
        self.notebook = self.interface.notebook
        self.add(self.notebook)
        # connect to the window's delete event to close on x click
        self.connect('destroy', gtk.main_quit)
        # set up indicator of changes, tabs are not disabled initially
        self._change_indicator = [False] * self.notebook.get_n_pages()
        self._reading_plates = True
        self.regex_model, self.wellsplate_model = gtk.ListStore(str, str, str, str), gtk.ListStore(str, str)

        # general window setup
        ##########################

        self._init_FILES()
        self._init_REGEX()

    ##########################
    # MARK functions general

    def _set_page_changed(self, changed):
        """
        Set page to changed and disable further stages; set page to unchanged and enable next page
        :param changed: toggle
        :return:
        """
        page = self.notebook.get_current_page()
        if changed:
            # disable later pages
            [page.set_sensitive(False) for page in
             self.notebook.get_children()[page + 1: self.notebook.get_n_pages()]]
            # set later pages to 'changed'
            self._change_indicator[page:] = [True] * (self.notebook.get_n_pages() - page)
        else:
            # enable the next page
            self.notebook.get_children()[page + 1].set_sensitive(True)
            # set earlier pages to 'unchanged'
            self._change_indicator[:page + 1] = [False] * (page + 1)

    def _get_page_changed(self):
        return self._change_indicator[self.notebook.get_current_page()]

    def _proceed(self, widget):
        # first integrate changes to the dataset
        if self._get_page_changed():
            self._prep_next(self.notebook.get_current_page())
        # then proceed
        self.notebook.next_page()
        self.log.debug('proceeded to page %d' % self.notebook.get_current_page())

    def _step_back(self, widget):
        self.notebook.prev_page()
        self._set_page_changed(False)
        self.log.debug('stepped back to page %d' % self.notebook.get_current_page())

    def _prep_next(self, page):
        """
        TODO page-specific handlers
        :param self:
        :param page:
        :return:
        """
        if page == 0:
            self._reset_REGEX()

        # make changes in dataset
        # self.data.prep_next(page)
        # transfer changes from dataset to GUI?

        self._set_page_changed(False)
        return

    ##########################
    # DONE page FILES

    def _init_FILES(self):
        # trace file types
        # create a TreeView model
        self.trace_model = gtk.ListStore(str, bool)
        [self.trace_model.append([filetype, False]) for filetype in gui.FILETYPES]

        # check ABI traces by default
        self.trace_model[0][1] = True
        self.data.filetypes.add('.ab1')

        self.interface.view_filetypes.set_model(self.trace_model)
        self.interface.view_filetypes.set_headers_visible(False)
        self.interface.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Filetype', cell_renderer=gtk.CellRendererText(), text=0))
        crt = gtk.CellRendererToggle()
        crt.connect('toggled', self._change_filetypes, self.data.filetypes, self.trace_model)
        self.interface.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))

        # trace file paths
        # create a TreeView model
        self.trace_model = gtk.ListStore(str)
        self.interface.view_trace_path.set_model(self.trace_model)
        self.interface.view_trace_path.set_headers_visible(False)
        self.interface.view_trace_path.append_column(
            gtk.TreeViewColumn(title='Paths', cell_renderer=gtk.CellRendererText(), text=0))
        trace_selection = self.interface.view_trace_path.get_selection()
        trace_selection.set_mode(gtk.SelectionMode.MULTIPLE)

        # connect buttons
        self.interface.add_trace_folder.connect('clicked', self._add_folder,
                                                self.data.filetypes, self.data.trace_paths, self.trace_model)
        self.interface.add_trace_manual.connect('clicked', self._add_manually,
                                                self.data.trace_paths, self.trace_model)
        self.interface.add_trace_whitelist.connect('clicked', self._add_manually,
                                                   self.data.trace_paths, self.trace_model, True)
        self.interface.remove_traces.connect('clicked', self._remove_path, trace_selection, self.data.trace_paths)
        self.interface.delete_all_traces.connect('clicked', self._clear_path, self.data.trace_paths, self.trace_model)

        # wellsplate paths
        # create a TreeView model
        self.csv_model = gtk.ListStore(str)
        self.interface.view_csv_path.set_model(self.csv_model)
        self.interface.view_csv_path.set_headers_visible(False)
        self.interface.view_csv_path.append_column(
            gtk.TreeViewColumn(title='Paths', cell_renderer=gtk.CellRendererText(), text=0))
        csv_selection = self.interface.view_csv_path.get_selection()
        csv_selection.set_mode(gtk.SelectionMode.MULTIPLE)

        # connect buttons
        self.interface.add_csv_folder.connect('clicked', self._add_folder,
                                              {'.csv'}, self.data.csv_paths, self.csv_model)
        self.interface.add_csv_manual.connect('clicked', self._add_manually, self.data.csv_paths, self.csv_model)
        self.interface.remove_csv.connect('clicked', self._remove_path,
                                          csv_selection, self.data.csv_paths)
        self.interface.delete_all_csv.connect('clicked', self._clear_path, self.data.csv_paths, self.csv_model)

        self.interface.files_next.connect('clicked', self._proceed)

    def _change_filetypes(self, widget, path, filetypes, model):
        # toggle the button
        model[path][1] = not model[path][1]
        # adapt the filetypes to read for the next folder
        if model[path][1]:
            filetypes.add(model[path][0])
        else:
            filetypes.remove(model[path][0])
        self.log.debug('selected ' + ' and '.join(filetypes))

    def _add_folder(self, widget, types, store, model):
        dialog = gtk.FileChooserDialog(title='select folder(s)', parent=self, select_multiple=True,
                                       action=gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                           gtk.STOCK_OPEN, gtk.ResponseType.OK)
        duplicates, new_paths = list(), list()
        response = dialog.run()
        if response == gtk.ResponseType.OK:
            for folder in dialog.get_filenames():
                try:
                    folder = Path(folder).resolve()
                    # get all the matching files
                    new_paths = [list(folder.glob('**/*%s' % ext)) for ext in types]
                    # flatten the LoL
                    new_paths = [str(item) for sublist in new_paths for item in sublist]
                    # save new entries in the dataset
                    store, new_paths, duplicates = self._add_new_entries(store, new_paths)
                    # append to the ListStore
                    [model.append([path]) for path in new_paths]
                except UnicodeDecodeError as ex:
                    self.log.info(ex)
            self._refresh_trace_number()
            self.log.info('found %d new paths in folder(s)' % len(new_paths))

        elif response == gtk.ResponseType.CANCEL:
            self.log.info('cancel')
        elif response == gtk.ResponseType.DELETE_EVENT:
            self.log.info('escape')

        dialog.destroy()

        if duplicates:
            self._show_message_dialog('Some files already selected', duplicates)

    def _add_manually(self, widget, store, model, whitelists=False):
        dialog = gtk.FileChooserDialog(title='select files', parent=self,
                                       action=gtk.FileChooserAction.OPEN, select_multiple=True)
        dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                           gtk.STOCK_OPEN, gtk.ResponseType.OK)
        not_found, duplicates = list(), list()
        response = dialog.run()
        if response == gtk.ResponseType.OK:
            # read the filenames
            new_paths = dialog.get_filenames()
            if whitelists:
                # read lines from the files
                extracted = list()
                for whitelist in new_paths:
                    try:
                        extracted += open(whitelist, 'r').read().strip().split('\n')
                    except UnicodeDecodeError as ex:
                        self.log.info(ex)

                # sort extracted into valid and invalid file paths
                new_paths = list()
                for string_path in extracted:
                    path = Path(string_path).resolve()
                    if path.is_file():
                        new_paths.append(str(path))
                    else:
                        not_found.append(string_path)

            # save new entries in the dataset
            store, new_paths, duplicates = self._add_new_entries(store, new_paths)
            # append to the ListStore
            [model.append([path]) for path in new_paths]
            self._refresh_trace_number()
            self.log.info('added %d paths' % len(new_paths))

        elif response == gtk.ResponseType.CANCEL:
            self.log.info('cancel')
        elif response == gtk.ResponseType.DELETE_EVENT:
            self.log.info('escape')

        dialog.destroy()

        if not_found:
            self._show_message_dialog('Some file paths were invalid', not_found)
        if duplicates:
            self._show_message_dialog('Some files already selected', duplicates)

    def _remove_path(self, widget, selection, store):
        model, iter = selection.get_selected_rows()
        [store.pop(item[0]) for item in reversed(iter)]
        model.clear()
        # re-append to the ListStore
        [model.append([path]) for path in store]
        self._set_page_changed(True)
        self._refresh_trace_number()
        self.log.info('removed %d rows' % len(iter))

    def _clear_path(self, widget, store, model):
        store.clear()
        model.clear()
        self._set_page_changed(True)
        self._refresh_trace_number()
        self.log.info('cleared paths')

    def _refresh_trace_number(self):
        """
        Adjusts the display of the number of selected trace files
        and toggles reading wellsplates or not.
        :return:
        """
        num_traces = len(self.trace_model)
        if num_traces > 0:
            self.interface.trace_number.set_label('%d files' % num_traces)
            self.interface.trace_number.set_visible(True)
        else:
            self.interface.trace_number.set_label('0 files')
            self.interface.trace_number.set_visible(False)

        # toggle _reading_plates
        self._reading_plates = len(self.csv_model) > 0
        self._set_page_changed(not self._reading_plates or self._get_page_changed())
        [self.interface.__getattribute__(name).set_sensitive(self._reading_plates) for name in ['plate_regex_label',
                                                                                                'plate_regex',
                                                                                                'wellsplate_regex_description',
                                                                                                'wellsplate_regex',
                                                                                                'wellsplate_buttons',
                                                                                                'wellsplate_regex_box']]
        # toggle radiobutton line back off
        if not self.interface.triple_regex_toggle.get_active():
            [self.interface.__getattribute__(name).set_sensitive(False) for name in
             ['plate_regex', 'plate_regex_label']]

    def _show_message_dialog(self, message, list_to_print=None):
        dialog = gtk.MessageDialog(transient_for=self, flags=0,
                                   buttons=gtk.ButtonsType.OK, message_type=gtk.MessageType.WARNING,
                                   text=message)
        # dialog.format_secondary_text('\n'.join(not_found))  # not copyable -> not user-friendly
        if list_to_print:
            txt_buf = gtk.TextBuffer()
            txt_buf.set_text('\n'.join(list_to_print))
            dialog.get_message_area().add(gtk.TextView().new_with_buffer(txt_buf))
        dialog.show_all()  # important
        dialog.run()
        dialog.destroy()
        self.log.debug('showed a message')

    def _add_new_entries(self, store, entries_to_add):
        """
        Adds new entries to a list and returns the modified list
        as well as lists of both newly added entries and observed duplicates
        :param store: the trace_paths list of a dataset object
        :param entries_to_add: obvs
        :return:
        """
        dups, news = list(), list()
        for entry in entries_to_add:
            if entry in store:
                dups.append(entry)
            else:
                store.append(entry)
                news.append(entry)
        if len(news) > 0:
            self._set_page_changed(True)
        return store, news, dups

    ##########################
    # TODO page RegEx

    def _init_REGEX(self):
        # connect buttons
        self.interface.single_regex_toggle.connect('toggled', self._regex_toggle)
        self.interface.triple_regex_toggle.connect('toggled', self._regex_toggle)
        self.interface.triple_regex_toggle.join_group(self.interface.single_regex_toggle)
        self.interface.single_regex_toggle.set_active(True)

        self.interface.regex_apply.connect('clicked', self._parse)
        self.interface.single_regex.connect('activate', self._parse)
        self.interface.well_regex.connect('activate', self._parse_single, self.interface.well_regex, 0)
        self.interface.gene_regex.connect('activate', self._parse_single, self.interface.gene_regex, -2)
        self.interface.plate_regex.connect('activate', self._parse_single, self.interface.plate_regex, 1)

        self.interface.wellsplate_regex.connect('activate', self._parse_single, self.interface.wellsplate_regex, 0)
        self.interface.wellsplate_apply.connect('clicked', self._parse_single, self.interface.wellsplate_regex, 0)

        # TODO on-the-fly RegEx help page
        # TODO save parser results in data

        # connect buttons
        self.interface.regex_next.connect('clicked', self._proceed)
        self.interface.regex_back.connect('clicked', self._step_back)
        self._reset_REGEX()
        self._refresh_trace_number()

    def _reset_REGEX(self):
        # path extraction: well/id, (plate), gene, file
        # create TreeView models and fill filename column
        if self._reading_plates:
            self.regex_model = gtk.ListStore(str, str, str, str)
            self.wellsplate_model = gtk.ListStore(str, str)
            # fill with initial data
            [self.regex_model.append([''] * 3 + [Path(path).name]) for path in self.data.trace_paths]
            [self.wellsplate_model.append([''] + [Path(path).name]) for path in self.data.csv_paths]
        else:
            self.regex_model = gtk.ListStore(str, str, str)
            [self.regex_model.append([''] * 2 + [Path(path).name]) for path in self.data.trace_paths]
        # self.interface.view_trace_regex.set_headers_visible(False)
        self.interface.view_trace_regex.set_model(self.regex_model)
        self.interface.view_csv_regex.set_model(self.wellsplate_model)
        # remove old columns:
        for treeview in [self.interface.view_trace_regex, self.interface.view_csv_regex]:
            [treeview.remove_colum(col) for col in treeview.get_columns()]

        # columns depend on _reading_plates
        self.interface.view_trace_regex.get_columns().clear()
        if self._reading_plates:
            for title, column in zip(['well', 'plate', 'gene', 'file'], list(range(4))):
                self.interface.view_trace_regex.append_column(
                    gtk.TreeViewColumn(title=title, cell_renderer=gtk.CellRendererText(), markup=column))
            # wellsplates:
            for title, column in zip(['plate ID', 'file'], list(range(2))):
                self.interface.view_csv_regex.append_column(
                    gtk.TreeViewColumn(title='title', cell_renderer=gtk.CellRendererText(), markup=column))
        else:
            for title, column in zip(['sample', 'gene', 'file'], list(range(3))):
                self.interface.view_trace_regex.append_column(
                    gtk.TreeViewColumn(title=title, cell_renderer=gtk.CellRendererText(), markup=column))

        for treeview in [self.interface.view_trace_regex, self.interface.view_csv_regex]:
            treeview.columns_autosize()
        for col_index in range(self.interface.view_trace_regex.get_n_columns()):
            self.interface.view_trace_regex.get_column(col_index).set_sort_column_id(col_index)

    def _regex_toggle(self, widget):
        # (In)activates the Entry fields and causes a parse event
        if widget.get_active():
            three = ['plate_regex',
                     'gene_regex',
                     'well_regex',
                     'plate_regex_label',
                     'gene_regex_label',
                     'well_regex_label']
            if widget == self.interface.single_regex_toggle:
                self.interface.single_regex.set_sensitive(True)
                [self.interface.__getattribute__(widget).set_sensitive(False) for widget in three]
            else:
                self.interface.single_regex.set_sensitive(False)
                [self.interface.__getattribute__(widget).set_sensitive(True) for widget in three]
            # adjust the sensitivity of the plate regex line
            if not self._reading_plates:
                self.interface.plate_regex.set_sensitive(False)
                self.interface.plate_regex_label.set_sensitive(False)
            self.log.debug('toggled radiobutton %s' % widget.get_name())
            # parse again
            self._parse(widget)

    def _parse(self, widget):
        if self.interface.single_regex_toggle.get_active():
            self.log.debug('parsing with single regex')
            # if parsing with only a single regex
            regex = re.compile(self.interface.single_regex.get_text())

            for row_index in range(len(self.regex_model)):
                file = self.regex_model[row_index][-1]
                try:
                    m = regex.search(file)
                    plate, gene, well = m.groups() if self._reading_plates else (None, *m.groups())
                    self.regex_model[row_index] = [well, plate, gene, file] if self._reading_plates else [well, gene,
                                                                                                          file]
                except ValueError as ve:
                    self.regex_model[row_index][-2] = '<span foreground="blue">groups?</span>'
                except AttributeError as ae:
                    self.regex_model[row_index][-2] = '<span foreground="red">no match</span>'
        else:
            self._parse_single(None, self.interface.well_regex, 0)
            self._parse_single(None, self.interface.gene_regex, -2)
            if self._reading_plates:
                self._parse_single(None, self.interface.plate_regex, 1)
        self.interface.view_trace_regex.columns_autosize()
        self.log.debug('parsing done')

    def _parse_single(self, widget, entry, col_index):
        self.log.debug('parsing from %s' % entry.get_name())
        regex = re.compile(entry.get_text())

        if widget in [self.interface.wellsplate_regex, self.interface.wellsplate_apply]:
            model = self.wellsplate_model
        else:
            model = self.regex_model

        for row_index in range(len(model)):
            file = model[row_index][-1]
            try:
                model[row_index][col_index] = regex.search(file).groups()[0]
            except ValueError as ve:
                model[row_index][col_index] = '<span foreground="blue">groups?</span>'
            except AttributeError as ae:
                model[row_index][col_index] = '<span foreground="red">no match</span>'


class dataset:
    def __init__(self):
        self.filetypes = set()
        self.trace_paths = []
        self.csv_paths = []
        self.path_container = []


def _init_log(**kwargs):
    """Initializes logging."""
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)

    if 'filename' in kwargs:
        # init verbose logging to file
        fh = logging.FileHandler(filename=kwargs['filename'], mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s\t%(name)s\t%(message)s',
                                          datefmt='%Y-%m-%d %H:%M:%S'))
        log.addHandler(fh)

    # init shortened console logging
    sh = logging.StreamHandler(sys.stdout)
    if __verbose__:
        sh.setLevel(logging.DEBUG)
    elif __info__:
        sh.setLevel(logging.INFO)
    else:
        sh.setLevel(logging.WARNING)
    sh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    log.addHandler(sh)


_init_log()  # filename='nope')
log = logging.getLogger(__name__)
log.debug('AB12PHYLO GUI version')

data = dataset()
win = gui(data)
win.show_all()
gtk.main()
