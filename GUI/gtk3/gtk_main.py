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
__verbose__, __info__ = 1,0

# set the icon theme
gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


# http://cdn.php-gtk.eu/cdn/farfuture/riUt0TzlozMVQuwGBNNJsaPujRQ4uIYXc8SWdgbgiYY/mtime:1368022411/sites/php-gtk.eu/files/gtk-php-get-child-widget-by-name.php__0.txt
# note get_name() vs gtk.Buildable.get_name(): https://stackoverflow.com/questions/3489520/python-gtk-widget-name
def get_descendant(widget, child_name, level, doPrint=False):
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
            return get_descendant(child, child_name, level + 1, doPrint)
    # It might have many children, so search them
    elif hasattr(widget, 'get_children') and callable(getattr(widget, 'get_children')) and child_name != "":
        children = widget.get_children()
        # For each child
        found = None
        for child in children:
            if child is not None:
                found = get_descendant(child, child_name, level + 1, doPrint)  # //search the child
                if found:
                    return found


class gui(gtk.Window):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'
    ELEMENTS = ['add_trace_folder', 'add_trace_manual', 'add_trace_whitelist', 'add_csv_folder', 'add_csv_manual',
                'view_filetypes', 'view_trace_path', 'trace_number', 'files_next', 'remove_path', 'delete_all']
    FILETYPES = ['.ab1', '.seq', '.fasta']

    def __init__(self, data):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(640, 480)
        self.set_size_request(640, 480)
        self.log.debug('GTK Window initialized')

        self.data = data

        # fetch the notebook from the .glade XML
        notebook = gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_object('notebook')
        self.add(notebook)
        self.log.debug('Fetched notebook')

        # connect to the window's delete event to close on x click
        self.connect('destroy', gtk.main_quit)
        self.ns = Namespace(**{widget: get_descendant(self, widget, 10, False) for widget in gui.ELEMENTS})
        self.log.debug('Fetched control elements')

        # DONE trace file types
        # create a TreeView model
        self.ns.file_type_model = gtk.ListStore(str, bool)
        [self.ns.file_type_model.append([filetype, False]) for filetype in gui.FILETYPES]

        # check ABI traces by default
        self.ns.file_type_model[0][1] = True
        self.data.trace_types.add('.ab1')

        self.ns.view_filetypes.set_model(self.ns.file_type_model)
        self.ns.view_filetypes.set_headers_visible(False)
        self.ns.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Filetype', cell_renderer=gtk.CellRendererText(), text=0))
        crt = gtk.CellRendererToggle()
        crt.connect('toggled', self._change_trace_types)
        self.ns.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))

        # TODO trace file paths
        # create a TreeView model
        self.ns.trace_path_model = gtk.ListStore(str)
        self.ns.view_trace_path.set_model(self.ns.trace_path_model)
        self.ns.view_trace_path.set_headers_visible(False)
        self.ns.view_trace_path.append_column(
            gtk.TreeViewColumn(title='Paths', cell_renderer=gtk.CellRendererText(), text=0))
        self.trace_selection = self.ns.view_trace_path.get_selection()
        self.trace_selection.set_mode(gtk.SelectionMode.MULTIPLE)

        # self.
        # self.ns.path_toolbar.set_icon_size(gtk.IconSize.SMALL_TOOLBAR)

        # connect manual file selection
        self.ns.add_trace_manual.connect('clicked', self._add_traces_manually)
        # connect whitelisted paths
        self.ns.add_trace_whitelist.connect('clicked', self._add_whitelist_traces)
        # connect removing some paths
        self.ns.remove_path.connect('clicked', self._remove_trace_paths)
        # connect clearing paths
        self.ns.delete_all.connect('clicked', self._clear_trace_paths)

    def _change_trace_types(self, widget, path):
        # toggle the button
        self.ns.file_type_model[path][1] = not self.ns.file_type_model[path][1]
        # adapt the filetypes to read for the next folder
        if self.ns.file_type_model[path][1]:
            self.data.trace_types.add(self.ns.file_type_model[path][0])
        else:
            self.data.trace_types.remove(self.ns.file_type_model[path][0])
        self.log.debug('selected ' + ' and '.join(self.data.trace_types))

    def _add_traces_manually(self, widget):
        trace_dialog = gtk.FileChooserDialog(title='Select trace files', parent=self,
                                             action=gtk.FileChooserAction.OPEN,
                                             select_multiple=True)
        trace_dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                                 gtk.STOCK_OPEN, gtk.ResponseType.OK)
        try:
            response = trace_dialog.run()
            if response == gtk.ResponseType.OK:
                self.log.debug('manual trace dialog returned ok')
                new_paths = trace_dialog.get_filenames()

                # save in the dataset
                self.data.trace_paths += new_paths
                # append to the ListStore
                [self.ns.trace_path_model.append([trace]) for trace in new_paths]
                self._refresh_trace_number()

            elif response == gtk.ResponseType.CANCEL:
                raise ValueError('cancel')
            elif response == gtk.ResponseType.DELETE_EVENT:
                raise ValueError('escape')

        except ValueError as ex:
            self.log.info(ex)
        trace_dialog.destroy()

    def _add_whitelist_traces(self, widget):
        trace_dialog = gtk.FileChooserDialog(title='Select whitelist file', parent=self,
                                             action=gtk.FileChooserAction.OPEN,
                                             select_multiple=True)
        trace_dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                                 gtk.STOCK_OPEN, gtk.ResponseType.OK)
        try:
            response = trace_dialog.run()
            if response == gtk.ResponseType.OK:
                self.log.debug('whitelist dialog returned ok')
                whitelists = trace_dialog.get_filenames()

                for whitelist in whitelists:
                    try:
                        new_paths = open(whitelist, 'r').read().strip().split('\n')
                        # save in the dataset
                        self.data.trace_paths += new_paths
                        # append to the ListStore
                        [self.ns.trace_path_model.append([trace]) for trace in new_paths]
                    except UnicodeDecodeError as ex:
                        self.log.info(ex)
                    trace_dialog.destroy()
                self._refresh_trace_number()

            elif response == gtk.ResponseType.CANCEL:
                raise ValueError('cancel')
            elif response == gtk.ResponseType.DELETE_EVENT:
                raise ValueError('escape')

        except ValueError as ex:
            self.log.info(ex)


    def _remove_trace_paths(self, widget):
        model, iter = self.trace_selection.get_selected_rows()
        [self.data.trace_paths.pop(item[0]) for item in reversed(iter)]
        model.clear()
        # re-append to the ListStore
        [self.ns.trace_path_model.append([trace]) for trace in self.data.trace_paths]
        self.log.debug('removed %d trace paths' % len(iter))
        self._refresh_trace_number()

    def _clear_trace_paths(self, widget):
        self.data.trace_paths = []
        self.ns.trace_path_model.clear()
        self._refresh_trace_number()


    def _refresh_trace_number(self):
        num_traces = len(self.data.trace_paths)
        if num_traces > 0:
            self.ns.trace_number.set_label('%d files' % num_traces)
            self.ns.trace_number.set_visible(True)
        else:
            self.ns.trace_number.set_label('0 files')
            self.ns.trace_number.set_visible(False)
        self.log.debug('%d trace files' % num_traces)


class dataset:
    def __init__(self):
        self.trace_types = set()
        self.trace_paths = []
        self.csv_paths = []


class driver:

    def __init__(self):
        self._init_log()  # filename='nope')
        log = logging.getLogger(__name__)
        log.debug('AB12PHYLO GUI version')

        data = dataset()
        win = gui(data)

        # # fetch the main window from the .glade XML
        # builder = gtk.Builder()
        # builder.add_from_file(str(TEMPLATE))
        # self.win = builder.get_object('win')
        # # self.win = gtk.Builder().new_from_file(str(TEMPLATE)).get_object('win')
        #
        # # plus = gtk.IconTheme().get_default().load_icon('list-add', 48, gtk.IconLookupFlags.NO_SVG)
        # # clicker = gtk.Button().new_from_icon_name('gtk-add', gtk.IconSize.BUTTON)  #.STOCK_GO_FORWARD, gtk.IconSize.BUTTON)
        # # # clicker.connect('clicked', self.on_button_clicked)
        # # builder.get_object('gridder').attach(clicker, 2, 2, 1, 1)
        #

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
