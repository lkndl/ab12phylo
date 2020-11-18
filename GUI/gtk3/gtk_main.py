# 2020 Leo Kaindl

import logging
import sys
import pickle
import threading
from argparse import Namespace
from pathlib import Path

import gi

from GUI.gtk3 import commons, files, regex, quality, align

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
        # self.set_size_request(1000, 800)

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
        iface.running = False
        self.project_path = None

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
        commons.bind_accelerator(self.accelerators, iface.open, '<Control>o', 'activate')
        iface.save.connect('activate', self.save)
        commons.bind_accelerator(self.accelerators, iface.save, '<Control>s', 'activate')
        iface.saveas.connect('activate', self.saveas)
        commons.bind_accelerator(self.accelerators, iface.saveas, '<Control><Shift>s', 'activate')

        self.data = project_dataset(iface.notebook.get_n_pages())
        self.log.debug('vars and dataset initialized')

        # initialize the notebook pages
        files.init(self)
        regex.init(self)
        quality.init(self)
        align.init(self)

        self.load('/home/quirin/PYTHON/AB12PHYLO/projects/stam.proj')

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

    def open(self, event):
        """
        Project file open dialog. Calls self.load internally.
        :param event: will be ignored, but must have an additional positional argument for use as callback
        :return:
        """
        dialog = Gtk.FileChooserDialog(title='open project', parent=None,
                                       action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            self.load(dialog.get_filename())
            dialog.destroy()
            return True
        dialog.destroy()
        return False

    def load(self, path):
        """
        Load a project.
        :param path: a string file path for the project file
        :return:
        """
        self.project_path = Path(path)
        self.log.debug('got dataset path %s' % self.project_path)
        # read in dataset
        with open(self.project_path, 'rb') as proj:
            self.data = pickle.load(proj)
        # set models again. do not delete next line, otherwise data refers to old one
        data, iface = self.data, self.interface
        for mo, tv in zip([data.trace_store, data.plate_store,
                           data.trace_store, data.plate_store, data.qal_model],
                          [iface.view_trace_path, iface.view_csv_path,
                           iface.view_trace_regex, iface.view_csv_regex, iface.view_qal]):
            tv.set_model(mo)
        iface.notebook.set_current_page(data.page)
        # set gene chooser + plot quality
        quality.reset(self)

    def save(self, event):
        if not self.project_path:
            self.saveas(None)
            return
        with open(self.project_path, 'wb') as proj:
            self.log.info('saving to %s' % self.project_path)
            try:
                self.data.page = self.interface.notebook.get_current_page()
                pickle.dump(self.data, proj)
                self.log.info('finished save')
            except pickle.PicklingError as pe:
                self.log.warning('saving failed')

    def saveas(self, event):
        dialog = Gtk.FileChooserDialog(title='save project', parent=None,
                                       action=Gtk.FileChooserAction.SAVE)
        if not self.project_path:
            dialog.set_current_name('untitled.proj')
        else:
            dialog.set_filename(str(self.project_path))
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            self.project_path = Path(dialog.get_filename())
            self.log.debug('got save path to %s' % self.project_path)
            self.save(None)
        dialog.destroy()


class project_dataset:
    def __init__(self, n_pages):
        self.trace_store = commons.picklable_liststore(str,  # path
                                                       str,  # filename
                                                       str,  # well/id
                                                       str,  # plate
                                                       str,  # gene
                                                       bool,  # reference
                                                       bool,  # reversed
                                                       str)  # color

        self.plate_store = commons.picklable_liststore(str,  # path
                                                       str,  # filename
                                                       str,  # plate ID
                                                       str)  # errors
        self.genes = set()  # used *before* seqdata exists
        self.csvs = dict()
        self.seqdata = dict()
        self.metadata = dict()
        self.seed = 0
        self.record_order = list()
        self.qal_model = commons.picklable_liststore(str,  # id
                                                     bool,  # has phreds
                                                     bool)  # low quality

        # set up indicator of changes, tabs are not disabled initially
        self.change_indicator = [False] * n_pages
        self.errors_indicator = [False] * n_pages
        self.page = 0

    def reset(self):
        for attr in [a for a in dir(self) if not callable(getattr(self, a))]:
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
