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
TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'
ELEMENTS = ['add_trace_folder', 'add_trace_manual', 'add_trace_whitelist', 'add_csv_folder', 'add_csv_manual',
            'view_filetypes', 'view_trace_path', 'trace_number', 'files_next', 'path_toolbar']
__verbose__, __info__ = 0, 1

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

    def __init__(self):
        log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(ICON))
        self.set_default_size(640, 480)
        self.set_size_request(640, 480)
        log.debug('GTK Window initialized')

        # fetch the notebook from the .glade XML
        notebook = gtk.Builder().new_from_file(str(TEMPLATE)).get_object('notebook')
        self.add(notebook)
        log.debug('Fetched notebook')

        # connect to the window's delete event to close on x click
        self.connect('destroy', gtk.main_quit)
        self.ns = Namespace(**{widget: get_descendant(self, widget, 10, False) for widget in ELEMENTS})
        log.debug('Fetched control elements')

        # create a TreeView model
        self.ns.filetypes = gtk.ListStore(str, bool)
        [self.ns.filetypes.append([filetype, False]) for filetype in ['.ab1', '.seq', '.fasta']]
        # check ABI traces by default
        self.ns.filetypes[0][1] = True
        self.ns.view_filetypes.set_model(self.ns.filetypes)
        self.ns.view_filetypes.set_headers_visible(False)
        self.ns.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Filetype', cell_renderer=gtk.CellRendererText(), text=0))
        crt = gtk.CellRendererToggle()
        crt.connect('toggled', self._on_filetype_toggled)
        self.ns.view_filetypes.append_column(
            gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))


        # TODO trace file paths
        # self.ns.path_toolbar.set_icon_size(gtk.IconSize.SMALL_TOOLBAR)

        # connect manual file selection
        self.ns.add_trace_manual.connect('clicked', self._add_traces_manually)


    def _on_filetype_toggled(self, widget, path):
        print(path)
        self.ns.filetypes[path][1] = not self.ns.filetypes[path][1]

    def _add_traces_manually(self, widget):
        print('yup')




class driver:

    def __init__(self):
        self._init_log()  # filename='nope')
        log = logging.getLogger(__name__)
        log.debug('AB12PHYLO GUI version')

        win = gui()

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

        self.csv_paths = []
        self.trace_paths = []
        self.trace_dirs = []
        self.trace_types = []

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
