import logging
import sys
from argparse import Namespace
from pathlib import Path

import gi

from GUI.gtk3 import files, regex, quality, commons

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib

BASE_DIR = Path(__file__).resolve().parents[2]
__verbose__, __info__ = 1, 0

# set the icon theme
Gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
Gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


class gui(Gtk.Window):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'

    def __init__(self):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        Gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(900, 600)
        self.set_size_request(640, 480)
        self.log.debug('GTK Window initialized')

        # fetch all named objects from the .glade XML
        self.interface = dict()
        for widget in Gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                self.interface[widget.get_name()] = widget
        self.interface = Namespace(**self.interface)
        self.log.debug('Fetched control elements')

        # get some colors
        sc = self.get_style_context()
        self.interface.BG = sc.get_background_color(Gtk.StateType.NORMAL)
        self.interface.FG = sc.get_color(Gtk.StateType.NORMAL)
        self.interface.BLUE = Gdk.RGBA(0.137255, 0.454902, 0.686275, 1)  # '#2374AF'
        self.interface.AQUA = Gdk.RGBA(0.180392, 0.701961, 0.596078, 1)  # '#2EB398'
        self.interface.RED = Gdk.RGBA(1, 0, 0, 1)

        # fetch the notebook
        self.notebook = self.interface.notebook
        self.add(self.notebook)
        # connect to the window's delete event to close on x click
        self.connect('destroy', Gtk.main_quit)

        # set up indicator of changes, tabs are not disabled initially
        self.interface.change_indicator = [False] * self.notebook.get_n_pages()
        self.interface.errors_indicator = [False] * self.notebook.get_n_pages()

        self.data = dataset()

        files.init(self)
        # regex.init(self)
        # quality.init(self)


class dataset:
    def __init__(self):
        self.filetypes = set()

        # paths and extracted data
        self.trace_store = Gtk.ListStore(str,  # path
                                         str,  # filename
                                         str,  # well/id
                                         str,  # plate
                                         str,  # gene
                                         bool,  # reference
                                         bool,  # reversed
                                         Gdk.RGBA)  # color

        self.plate_store = Gtk.ListStore(str,  # path
                                         str,  # filename
                                         str)  # plate ID

        # self.trace_model = Gtk.ListStore(str, bool, Gdk.RGBA)
        # self.csv_model = Gtk.ListStore(str, bool, Gdk.RGBA)
        # self.rx_model = Gtk.ListStore(str, str, str, str, bool, str, str)
        # self.wp_model = Gtk.ListStore(str, str, str)
        #
        # self.q_model = Gtk.ListStore(str, str)
        self.genes = list()
        self.csvs = dict()
        self.seqdata = dict()


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
log.info('AB12PHYLO GUI version')

win = gui()
win.show_all()
Gtk.main()