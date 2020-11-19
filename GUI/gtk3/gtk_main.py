# 2020 Leo Kaindl

import logging
import sys
import pickle
import threading
import warnings
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

    KXLIN = [('C', (0.46, 1, 0.44, 1)),
             ('G', (0.16, 0.44, 0.8, 1)),
             ('T', (1, 0.47, 0.66, 1)),
             ('A', (0.92, 1, 0.4, 1)),
             ('N', (0.84, 0.84, 0.84, 0.6)),
             ('-', (1, 1, 1, 0)),
             ('[ ]', (1, 1, 1, 0)),
             ('sep', (1, 1, 1, 0)),
             ('unknown', (1, 0, 0, 1))]

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
        self.set_hide_titlebar_when_maximized(True)
        self.log.debug('GTK Window initialized')

        # set up an empty thread so Gtk can get used to it
        iface.thread = threading.Thread()
        iface.running = False
        self.project_path = None
        self.wd = None

        # get some CSS styling
        mod = b'lighter', b'darker'
        mod2 = 200  # per default, make treeview text color darker
        if 'dark' in Gtk.Settings.get_default().get_property('gtk-theme-name'):
            mod = mod[::-1]
            mod2 = 255

        css_provider = Gtk.CssProvider()
        css_provider.load_from_data(b'''
        .codeview text { background-color: %s(@bg_color); color: %s(@fg_color); }''' % mod)
        # treeview {background-color: darker(@bg_color);}
        Gtk.StyleContext.add_provider_for_screen(Gdk.Screen.get_default(), css_provider,
                                                 Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION)
        # get some colors
        iface.RED = '#FF0000'
        iface.BLUE = '#2374AF'
        iface.GREEN = '#23AF46'
        iface.AQUA = '#2EB398'
        sc = self.get_style_context()
        iface.FG = '#' + ''.join([(hex(min(int(c * mod2), 255))[2:]).upper()
                                  for c in list(sc.get_color(Gtk.StateFlags.ACTIVE))[:-1]])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            iface.BG = sc.get_background_color(Gtk.StateType.NORMAL)
        # replace white as the non-color with the background color. lighten it a bit to match better
        iface.colors = gui.KXLIN[0:5] \
                       + list(zip([k[0] for k in gui.KXLIN[5:-1]],
                                  [tuple(round(min(1, c * 1.2), 2) for c in iface.BG)] * 3)) \
                       + gui.KXLIN[-1:]

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

        self.data = project_dataset()
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
            if response != Gtk.ResponseType.OK:
                return
        LOG.debug('new project')
        self.data.new_project()
        files.refresh_files(self)
        self.interface.notebook.set_current_page(0)

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
        Load a project from a file, overwriting previous data in-place.
        :param path: a string file path to the project file to load
        :return:
        """
        self.project_path = Path(path)
        self.log.debug('got dataset path %s' % self.project_path)
        # read in dataset
        with open(self.project_path, 'rb') as proj:
            new_data = pickle.load(proj)
        # overwrite content in old dataset in-place rather than re-pointing everything
        self.data.overwrite(new_data)
        self.interface.notebook.set_current_page(self.data.page)
        # set gene chooser + plot quality
        quality.reset(self)
        self.wd = self.project_path.parent / self.project_path.stem
        Path.mkdir(self.wd, exist_ok=True)

    def save(self, event):
        """
        Save project_dataset to file directly, unless it hasn't previously been saved.
        :param event: To ignore, GTK+ callback requires positional argument.
        :return:
        """
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
        self.wd = self.project_path.parent / self.project_path.stem
        Path.mkdir(self.wd, exist_ok=True)

    def saveas(self, event):
        """
        Save with a file dialog. If previously saved, suggest old filename.
        :param event: To ignore, GTK+ callback requires positional argument.
        :return:
        """
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
    def __init__(self):
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
        self.change_indicator = [False] * 20
        self.errors_indicator = [False] * 20
        self.page = 0
        self.width = 0

    def new_project(self):
        self.overwrite(project_dataset())

    def overwrite(self, new_dataset):
        for attr in [a for a in dir(self) if not callable(getattr(self, a)) and not a.startswith('__')]:
            old = self.__getattribute__(attr)
            if type(old) == commons.picklable_liststore:
                old.clear()
                [old.append(row[:]) for row in new_dataset.__getattribute__(attr)]
            elif type(old) == dict:
                old.clear()
                old.update(new_dataset.__getattribute__(attr))
            else:
                try:
                    self.__setattr__(attr, new_dataset.__getattribute__(attr))
                except (AttributeError, TypeError) as ex:
                    LOG.exception(ex)


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
LOG.info('AB12PHYLO GUI version')

win = gui()
win.show_all()
Gtk.main()
