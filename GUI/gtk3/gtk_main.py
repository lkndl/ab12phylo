import logging
import sys
from argparse import Namespace
from pathlib import Path

import gi

from GUI.gtk3 import files, regex

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk

BASE_DIR = Path(__file__).resolve().parents[2]
__verbose__, __info__ = 1, 0

# set the icon theme
gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


class gui(gtk.Window):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'

    def __init__(self):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(900, 600)
        self.set_size_request(640, 480)
        self.log.debug('GTK Window initialized')

        # fetch all named objects from the .glade XML
        self.interface = dict()
        for widget in gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                self.interface[widget.get_name()] = widget
        self.interface = Namespace(**self.interface)
        self.log.debug('Fetched control elements')

        # fetch the notebook
        self.notebook = self.interface.notebook
        self.add(self.notebook)
        # connect to the window's delete event to close on x click
        self.connect('destroy', gtk.main_quit)
        # set up indicator of changes, tabs are not disabled initially

        self.interface.change_indicator = [False] * self.notebook.get_n_pages()
        self.interface.plates, self.interface.search_rev = True, False

        self.data = dataset()

        files.init(self.data, self.interface)
        regex.init(self.data, self.interface)


class dataset:
    def __init__(self):
        self.filetypes = set()
        self.trace_model, self.csv_model = gtk.ListStore(str), gtk.ListStore(str)
        self.rx_model, self.wp_model = gtk.ListStore(str, str, str, str, str), gtk.ListStore(str, str)



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

win = gui()
win.show_all()
gtk.main()
