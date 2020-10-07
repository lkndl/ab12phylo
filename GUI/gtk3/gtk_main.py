import sys

import gi
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


# http://cdn.php-gtk.eu/cdn/farfuture/riUt0TzlozMVQuwGBNNJsaPujRQ4uIYXc8SWdgbgiYY/mtime:1368022411/sites/php-gtk.eu/files/gtk-php-get-child-widget-by-name.php__0.txt
# note get_name() vs gtk.Buildable.get_name(): https://stackoverflow.com/questions/3489520/python-gtk-widget-name
def _get_descendant(widget, child_name, level, doPrint=False):
    if widget is not None:
        if doPrint:
            try:
                print("-" * level + gtk.Buildable.get_name(widget) + " :: " + widget.get_name())
            except TypeError:
                print("-" * level + "<no_name>" + " :: " + widget.get_name())
    else:
        if doPrint:
            print("-" * level + "None")
        return None
    # If it is what we are looking for
    if gtk.Buildable.get_name(widget) == child_name:  # not widget.get_name() !
        return widget
    # If this widget has one child only search its child
    if hasattr(widget, 'get_child') and callable(getattr(widget, 'get_child')) and child_name != "":
        child = widget.get_child()
        if child is not None:
            return _get_descendant(child, child_name, level + 1, doPrint)
    # It might have many children, so search them
    elif hasattr(widget, 'get_children') and callable(getattr(widget, 'get_children')) and child_name != "":
        children = widget.get_children()
        # For each child
        found = None
        for child in children:
            if child is not None:
                found = _get_descendant(child, child_name, level + 1, doPrint)  # //search the child
                if found:
                    return found


def _add_new_entries(store, entries_to_add):
    dups, news = list(), list()
    for entry in entries_to_add:
        if entry in store:
            dups.append(entry)
        else:
            store.append(entry)
            news.append(entry)
    return store, news, dups


class gui(gtk.Window):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'
    ELEMENTS = ['add_trace_folder', 'add_trace_manual', 'add_trace_whitelist', 'add_csv_folder', 'add_csv_manual',
                'view_filetypes', 'view_trace_path', 'view_csv_path', 'trace_number', 'files_next',
                'remove_path', 'delete_all', 'remove_csv', 'delete_all_csv']
    FILETYPES = ['.ab1', '.seq', '.fasta']

    def __init__(self, data):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(640, 480)
        self.set_size_request(640, 480)
        self.log.debug('GTK Window initialized')

        # self.data = data

        # fetch the notebook from the .glade XML
        self.notebook = gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_object('notebook')
        self.add(self.notebook)
        self.log.debug('Fetched notebook')

        # connect to the window's delete event to close on x click
        self.connect('destroy', gtk.main_quit)
        self.interface = Namespace(**{widget: _get_descendant(self, widget, 10, False) for widget in gui.ELEMENTS})
        self.log.debug('Fetched control elements')

        ##########################
        # page FILES
        # trace file types
        # create a TreeView model
        self.trace_model = gtk.ListStore(str, bool)
        [self.trace_model.append([filetype, False]) for filetype in gui.FILETYPES]

        # check ABI traces by default
        self.trace_model[0][1] = True
        data.filetypes.add('.ab1')

        self.interface.view_filetypes.set_model(self.trace_model)
        self.interface.view_filetypes.set_headers_visible(False)
        self.interface.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Filetype', cell_renderer=gtk.CellRendererText(), text=0))
        crt = gtk.CellRendererToggle()
        crt.connect('toggled', self._change_filetypes, data.filetypes)
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
                                                data.filetypes, data.trace_paths, self.trace_model)
        self.interface.add_trace_manual.connect('clicked', self._add_manually,
                                                data.trace_paths, self.trace_model)
        self.interface.add_trace_whitelist.connect('clicked', self._add_manually,
                                                   data.trace_paths, self.trace_model, True)
        self.interface.remove_path.connect('clicked', self._remove_path, trace_selection, data.trace_paths)
        self.interface.delete_all.connect('clicked', self._clear_path, data.trace_paths, self.trace_model)

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
                                              {'.csv'}, data.csv_paths, self.csv_model)
        self.interface.add_csv_manual.connect('clicked', self._add_manually, data.csv_paths, self.csv_model)
        self.interface.remove_csv.connect('clicked', self._remove_path,
                                          csv_selection, data.csv_paths)
        self.interface.delete_all_csv.connect('clicked', self._clear_path, data.csv_paths, self.csv_model)

        self.interface.files_next.connect('clicked', self._switch_to_page, 1)

        # DONE page Files
        ##########################
        # TODO page RegEx

    def _change_filetypes(self, widget, path, filetypes):
        # toggle the button
        self.trace_model[path][1] = not self.trace_model[path][1]
        # adapt the filetypes to read for the next folder
        if self.trace_model[path][1]:
            filetypes.add(self.trace_model[path][0])
        else:
            filetypes.remove(self.trace_model[path][0])
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
                    store, new_paths, duplicates = _add_new_entries(store, new_paths)
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
            store, new_paths, duplicates = _add_new_entries(store, new_paths)
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

    def _remove_path(self, widget, select, store):
        model, iter = select.get_selected_rows()
        [store.pop(item[0]) for item in reversed(iter)]
        model.clear()
        # re-append to the ListStore
        [model.append([path]) for path in store]
        self._refresh_trace_number()

    def _clear_path(self, widget, store, model):
        store = list()
        model.clear()
        self._refresh_trace_number()
        self.log.info('cleared paths')

    def _refresh_trace_number(self):
        num_traces = len(self.trace_model)
        if num_traces > 0:
            self.interface.trace_number.set_label('%d files' % num_traces)
            self.interface.trace_number.set_visible(True)
        else:
            self.interface.trace_number.set_label('0 files')
            self.interface.trace_number.set_visible(False)

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

    def _switch_to_page(self, widget, to):
        self.notebook.set_current_page(to)
        self.log.info('moved to page %d' % to)

    # DONE page Files
    ##########################
    # TODO page RegEx


class dataset:
    def __init__(self):
        self.filetypes = set()
        self.trace_paths = []
        self.csv_paths = []


class driver:

    def __init__(self):
        self._init_log()  # filename='nope')
        log = logging.getLogger(__name__)
        log.debug('AB12PHYLO GUI version')

        data = dataset()
        win = gui(data)
        win.show_all()
        gtk.main()

    def _init_log(self, **kwargs):
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


driver()
