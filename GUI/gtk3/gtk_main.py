# 2020 Leo Kaindl

import logging
import sys
import pickle
import threading
from argparse import Namespace
from pathlib import Path

import gi

from GUI.gtk3 import files, regex, quality, commons

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib, GObject

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
__verbose__, __info__ = 1, 0

# set the icon theme
Gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
Gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


class gui(Gtk.ApplicationWindow):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'

    def __init__(self):
        self.log = logging.getLogger(__name__)
        # super(gui, self).__init__()
        Gtk.Window.__init__(self, title='AB12PHYLO')
        self.set_icon_from_file(str(gui.ICON))
        self.set_default_size(900, 600)
        self.set_size_request(640, 480)

        # fetch all named objects from the .glade XML
        iface = dict()
        for widget in Gtk.Builder().new_from_file(str(gui.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                iface[widget.get_name()] = widget
        self.interface = Namespace(**iface)
        iface = self.interface

        # populate the window
        self.add(iface.toplayer)
        self.log.debug('GTK Window initialized')

        # set up an empty thread so Gtk can get used to it
        iface.thread = threading.Thread()
        self.cwd = Path.cwd()
        self.proj = None

        # set up indicator of changes, tabs are not disabled initially
        iface.change_indicator = [False] * iface.notebook.get_n_pages()
        iface.errors_indicator = [False] * iface.notebook.get_n_pages()

        # get some colors
        iface.RED = '#FF0000'
        iface.BLUE = '#2374AF'
        iface.GREEN = '#23AF46'
        iface.AQUA = '#2EB398'
        sc = self.get_style_context()
        iface.FG = '#' + ''.join([(hex(min(int(c * 256), 255))[2:]).upper()
                                  for c in list(sc.get_color(Gtk.StateType.NORMAL))[:-1]])
        # iface.BG = sc.get_background_color(Gtk.StateType.NORMAL)

        # prepare shortcuts / accelerators
        self.accelerators = Gtk.AccelGroup()
        self.add_accel_group(self.accelerators)

        # connect to the window's delete event to close on x click
        self.connect('destroy', Gtk.main_quit)
        # connect menu events and shortcuts
        iface.quit.connect('activate', Gtk.main_quit)
        commons.bind_accelerator(self.accelerators, iface.quit, '<Control>q', 'activate')
        iface.new.connect('activate', self.new)
        commons.bind_accelerator(self.accelerators, iface.new, '<Control>n', 'activate')
        iface.open.connect('activate', self.open)
        commons.bind_accelerator(self.accelerators, iface.open, '<Control>n', 'activate')
        iface.save.connect('activate', self.save)
        commons.bind_accelerator(self.accelerators, iface.save, '<Control>n', 'activate')
        iface.saveas.connect('activate', self.saveas)
        commons.bind_accelerator(self.accelerators, iface.saveas, '<Control>n', 'activate')

        self.data = dataset()
        self.log.debug('vars and dataset initialized')

        # initialize the notebook pages
        files.init(self)
        regex.init(self)
        quality.init(self)

    def new(self, confirm=True):
        """
        Create a new project by resetting the dataset and some interface aspects.
        :return:
        """
        if confirm:
            message = 'Create new project, discard unsaved changes?'
            dialog = Gtk.MessageDialog(transient_for=None, flags=0, text=message,
                                       buttons=Gtk.ButtonsType.OK_CANCEL,
                                       message_type=Gtk.MessageType.QUESTION)
            dialog.show_all()  # important
            response = dialog.run()
            dialog.destroy()

            if response == Gtk.ResponseType.OK:
                LOG.debug('new project')
                self.data.reset()
                files.refresh_files(self)
        else:
            self.data.reset()
            files.refresh_files(self)

    def open(self):
        pass
    def save(self):
        if not self.proj:
            self.saveas()
        with open(self.proj) as proj:
            pickle.dump((self.data, self.interface.toplayer), proj)
    def saveas(self):
        dialog = Gtk.FileChooserDialog(title='save project', parent=None, select_multiple=True,
                                       action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        pass


class dataset:
    def __init__(self):
        self.trace_store = Gtk.ListStore(str,  # path
                                         str,  # filename
                                         str,  # well/id
                                         str,  # plate
                                         str,  # gene
                                         bool,  # reference
                                         bool,  # reversed
                                         str)  # color

        self.plate_store = Gtk.ListStore(str,  # path
                                         str,  # filename
                                         str,  # plate ID
                                         str)  # errors
        self.genes = set()  # used *before* seqdata exists
        self.csvs = dict()
        self.seqdata = dict()
        self.metadata = dict()
        self.seed = 0
        self.record_order = list()
        self.qal_model = Gtk.ListStore(str,  # id
                                       bool,  # has phreds
                                       bool)  # low quality

    def reset(self):
        for attr in [a for a in dir(self)
                     if not a.startswith('__') and not callable(getattr(self, a))]:
            try:
                self.__getattribute__(attr).clear()
            except AttributeError:
                self.__setattr__(attr, 0)


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
