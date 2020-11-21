# 2020 Leo Kaindl


__author__ = 'Leo Kaindl'
__email__ = 'leo.kaindl@tum.de'
__version__ = '0.3a.01'
__date__ = '20 November 2020'
__license__ = 'MIT'
__status__ = 'Alpha'

import logging
import sys
import pickle
import threading
import warnings
import shutil
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


class app(Gtk.Application):
    TEMPLATE = BASE_DIR / 'GUI' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'GUI' / 'files' / 'favi.ico'

    def do_activate(self):
        self.add_window(self.win)
        self.win.show_all()

    def __init__(self):
        Gtk.Application.__init__(self)

        # fetch all named objects from the .glade XML
        iface = dict()
        for widget in Gtk.Builder().new_from_file(str(app.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                iface[widget.get_name()] = widget
        iface = Namespace(**iface)
        self.iface = iface
        self.win = iface.win
        self.win.set_icon_from_file(str(app.ICON))

        tbar = Gtk.HeaderBar()
        tbar.set_has_subtitle(False)
        tbar.set_decoration_layout('menu:minimize,maximize,close')
        tbar.set_show_close_button(True)
        tbar.pack_start(iface.menu_bar)
        self.win.set_titlebar(tbar)
        self.win.set_hide_titlebar_when_maximized(True)
        LOG.debug('GTK Window initialized')

        # set up an empty thread so Gtk can get used to it
        iface.thread = threading.Thread()
        iface.running = False
        iface.frac = 0
        iface.txt = ''
        self.project_path = None
        self.wd = Path.cwd() / 'untitled'
        Path.mkdir(self.wd, exist_ok=True)

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
        sc = self.win.get_style_context()
        iface.FG = '#' + ''.join([(hex(min(int(c * mod2), 255))[2:]).upper()
                                  for c in list(sc.get_color(Gtk.StateFlags.ACTIVE))[:-1]])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            iface.BG = sc.get_background_color(Gtk.StateType.NORMAL)

        # prepare shortcuts / accelerators
        self.accelerators = Gtk.AccelGroup()
        self.win.add_accel_group(self.accelerators)

        # connect to the window's delete event to close on x click
        self.win.connect('destroy', Gtk.main_quit)
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

        # connect buttons
        iface.next.connect('clicked', commons.proceed, self)
        commons.bind_accelerator(self.accelerators, iface.next, '<Alt>Right')
        iface.back.connect('clicked', commons.step_back, self)
        commons.bind_accelerator(self.accelerators, iface.back, '<Alt>Left')
        iface.refresh.connect('clicked', commons.re_run, self)
        commons.bind_accelerator(self.accelerators, iface.back, '<Control>r')
        # connect gene switcher
        iface.gene_handler = iface.gene_roll.connect('changed', commons.select, self)
        # any page change
        iface.notebook.connect_after('switch-page', commons.refresh, self)

        iface.dismiss.connect('clicked', lambda *args: iface.revealer.set_reveal_child(False))
        commons.bind_accelerator(self.accelerators, iface.dismiss, 'Escape')
        commons.bind_accelerator(self.accelerators, iface.dismiss, 'Return')

        self.data = project_dataset()
        LOG.debug('interface and dataset initialized')

        # initialize the notebook pages
        files.init(self)
        regex.init(self)
        quality.init(self)
        align.init(self)

        # self.load('/home/quirin/PYTHON/AB12PHYLO/projects/stam.proj')

    def new(self, confirm=True):
        """
        Create a new project by resetting the dataset and some interface aspects.
        :return:
        """
        if confirm:
            message = 'New project, discard changes?'
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
        self.wd = Path.cwd() / 'untitled'
        self.project_path = None
        self.iface.notebook.set_current_page(0)
        self.win.set_title('AB12PHYLO [untitled]')

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
        LOG.debug('got dataset path %s' % self.project_path)
        # read in dataset
        with open(self.project_path, 'rb') as proj:
            new_data = pickle.load(proj)
        # overwrite content in old dataset in-place rather than re-pointing everything
        self.data.overwrite(new_data)
        self.iface.notebook.set_current_page(self.data.page)
        self.wd = self.project_path.parent / self.project_path.stem
        Path.mkdir(self.wd, exist_ok=True)
        self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)

    def save(self, event, copy_from=None):
        """
        Save project_dataset to file directly, unless it hasn't previously been saved.
        :param event: To ignore, GTK+ callback requires positional argument.
        :return:
        """
        if not self.project_path:
            self.saveas(None)
            return
        with open(self.project_path, 'wb') as proj:
            LOG.info('saving to %s' % self.project_path)
            try:
                self.data.page = self.iface.notebook.get_current_page()
                pickle.dump(self.data, proj)
                LOG.info('finished save')
            except pickle.PicklingError as pe:
                LOG.warning('saving failed')

        new_wd = self.project_path.parent / self.project_path.stem
        if new_wd != self.wd:
            # in fact saving inside a saveas. -> copy from previous dir
            # clear new location, would be overwritten anyway;
            # doing it this roundabout way to keep compatibility with python3.6
            try:
                shutil.rmtree(path=new_wd)
            except FileNotFoundError:
                pass
            shutil.copytree(src=self.wd, dst=new_wd)  # not dirs_exist_ok=True

            # delete old data if it was in the prelim directory
            if self.wd == Path.cwd() / 'untitled':
                shutil.rmtree(path=self.wd)

        self.wd = self.project_path.parent / self.project_path.stem
        self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)

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
            LOG.debug('got save path to %s' % self.project_path)
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
        self.height = 0
        self.msa_hash = ''

    def agene(self):
        gene = self.genes.pop()
        self.genes.add(gene)
        return gene

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
                    LOG.error(ex)


def _init_log(**kwargs):
    """Initializes logging."""
    log = logging.getLogger()
    log.setLevel(logging.INFO)

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

app = app()
exit_status = app.run(sys.argv)
sys.exit(exit_status)
