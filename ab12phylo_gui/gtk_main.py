#!/usr/bin/python3
# 2021 Leo Kaindl

import logging
import pickle
import shutil
import sys
import threading
import warnings
from argparse import Namespace
from pathlib import Path

import gi

from ab12phylo.__init__ import __version__, __author__, __email__

from ab12phylo_gui.static import BASE_DIR
from ab12phylo_gui import gtk_proj, shared, gtk_io, gtk_rgx, \
    gtk_qal, gtk_msa, gtk_gbl, gtk_blast, gtk_ml, gtk_tree

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, Gio, GdkPixbuf

LOG = logging.getLogger(__name__)
__verbose__, __info__ = 1, 0

# set the icon theme
Gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Maia')
Gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-sea')


class ab12phylo_app(Gtk.Application):
    TEMPLATE = BASE_DIR / 'ab12phylo_gui' / 'files' / 'gui.glade'
    ICON = BASE_DIR / 'ab12phylo_gui' / 'files' / 'favi.ico'

    def do_activate(self):
        self.add_window(self.win)
        self.win.show_all()

    def do_shutdown(self):
        # delete data if it was in the prelim directory
        if self.wd == Path('untitled'):  # Path.cwd() / 'untitled':
            LOG.info('shutdown: delete prelim data')
            shutil.rmtree(path=self.wd)
        Gtk.Application.do_shutdown(self)

    def do_startup(self):
        Gtk.Application.do_startup(self)
        # connect menu actions and shortcuts
        for act, acc in zip(['new', 'open', 'save', 'save_as', 'help', 'about', 'on_quit'],
                            ['n', 'o', 's', '<Shift>s', 'h', '<Shift>h', 'q']):
            action = Gio.SimpleAction.new(act)
            action.connect('activate', self.__getattribute__(act))  # can pass parameters here
            self.add_action(action)
            self.set_accels_for_action(detailed_action_name='app.%s' % act, accels=['<Control>%s' % acc])
            # self.add_accelerator(accelerator='<Control>h', action_name='app.help', parameter=None)

    def __init__(self):
        Gtk.Application.__init__(self, application_id='org.ab12phylo', flags=Gio.ApplicationFlags.FLAGS_NONE)

        # fetch all named objects from the .glade XML
        iface = dict()
        for widget in Gtk.Builder().new_from_file(str(ab12phylo_app.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                iface[widget.get_name()] = widget
        iface = Namespace(**iface)
        self.iface = iface
        self.win = iface.win
        self.win.set_icon_from_file(str(ab12phylo_app.ICON))

        self.win.set_titlebar(iface.tbar)
        self.win.set_hide_titlebar_when_maximized(True)
        self.win.set_title('AB12PHYLO [untitled]')
        LOG.debug('GTK Window initialized')

        self.log = logging.getLogger()
        self._init_log()

        # set up an empty thread so Gtk can get used to it
        iface.thread = threading.Thread()
        iface.running = False
        iface.i = 0
        iface.k = 1
        iface.frac = 0
        iface.text = 'idle'
        iface.run_after = []

        # whether to draw raster (fast) or vector images is a button state

        # set up preliminary working directory
        self.project_path = None
        # self.wd = Path.cwd() / 'untitled'
        self.wd = Path('untitled')
        Path.mkdir(self.wd, exist_ok=True)

        # get some CSS styling
        mod = b'lighter', b'darker'
        mod2 = 200  # per default, make treeview text color darker
        theme = Gtk.Settings.get_default().get_property('gtk-theme-name')
        if 'dark' in theme:
            mod = mod[::-1]
            mod2 = 255

        css_provider = Gtk.CssProvider()
        css = b'''
        .seqid { font-size: xx-small; }
        separator.wide { min-width: 18px }
        progressbar trough progress { 
            min-height: 6px; 
            border-radius: 1px; 
            background-color: @success_color; }
        progressbar trough { min-height: 6px; }
        .hint { opacity: .6; } 
        #prog_label { font-size: x-small }
        button:active { background-color: #17f }
        notebook header { background-color: transparent; }
        '''
        if 'atcha' in theme:
            css += b'''
            .codeview text { background-color: %s(@bg_color); color: %s(@fg_color); }
            window { background-color: mix(@bg_color, rgba(127,127,127,0), 0.05) }
            ''' % mod
        else:
            css += b'''.codeview text { background-color: white; color: black; }'''
        css_provider.load_from_data(css)
        # CSS also supports key bindings
        # treeview {background-color: darker(@bg_color);} separator has margin-left
        Gtk.StyleContext.add_provider_for_screen(
            Gdk.Screen.get_default(), css_provider,
            Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION)
        # get some colors
        iface.RED = '#FF0000'
        iface.BLUE = '#2374AF'
        iface.GREEN = '#23AF46'
        iface.AQUA = '#2EB398'
        iface.PURPLE = '#8A2BE2'
        iface.tempspace = Namespace()
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
        self.win.connect('destroy', lambda *args: self.quit())

        # connect buttons
        iface.next.connect('clicked', shared.proceed, self)
        shared.bind_accelerator(self.accelerators, iface.next, '<Alt>Right')
        iface.back.connect('clicked', shared.step_back, self)
        shared.bind_accelerator(self.accelerators, iface.back, '<Alt>Left')
        iface.refresh.connect('clicked', shared.re_run, self)
        shared.bind_accelerator(self.accelerators, iface.refresh, 'Return')
        # connect gene switcher
        iface.gene_handler = iface.gene_roll.connect(
            'changed', shared.select_gene_and_redo, self)
        # any page change
        iface.notebook.connect_after('switch-page', shared.refresh, self)

        iface.dismiss.connect(
            'clicked', lambda *args: iface.revealer.set_reveal_child(False))
        shared.bind_accelerator(self.accelerators, iface.dismiss, 'Escape')
        shared.bind_accelerator(self.accelerators, iface.dismiss, 'Return')

        self.data = gtk_proj.project_dataset()
        LOG.debug('interface and dataset initialized')

        # initialize the notebook pages
        gtk_io.init(self)
        gtk_rgx.init(self)
        gtk_qal.init(self)
        gtk_msa.init(self)
        gtk_gbl.init(self)
        gtk_blast.init(self)
        gtk_ml.init(self)
        gtk_tree.init(self)

        self.load('ttest.proj')

    def new(self, action, confirm=True, *args):
        """
        Create a new project by resetting the dataset and some interface aspects.
        :return:
        """
        if confirm or confirm is None:
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
        # self.wd = Path.cwd() / 'untitled'
        self.wd = Path('untitled')
        self.project_path = None
        self.win.set_title('AB12PHYLO [untitled]')
        self.iface.notebook.set_current_page(0)

    def open(self, *args):
        """
        Project file open dialog. Calls self.load internally.
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
        project_path = Path(path)
        LOG.debug('got dataset path %s' % project_path)
        # read in dataset
        try:
            with open(project_path, 'rb') as proj:
                new_data = pickle.load(proj)
            # overwrite content in old dataset in-place rather than re-pointing everything
            self.data.overwrite(new_data)
            try:
                shared.init_gene_roll(self)
            except ValueError as ve:
                LOG.error(ve)
            self.project_path = project_path
            self.wd = self.project_path.parent / self.project_path.stem
            Path.mkdir(self.wd, exist_ok=True)
            self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)
            self.iface.notebook.set_current_page(self.data.page)

            for module in [gtk_proj, shared, gtk_io, gtk_rgx, gtk_qal,
                           gtk_msa, gtk_gbl, gtk_blast, gtk_ml, gtk_tree]:
                try:
                    module.reload_ui_state(self)
                except AttributeError:
                    pass  # TODO not elegant
        except Exception as e:
            LOG.exception(e)
            shared.show_notification(self, 'Project could not be loaded')

    def save(self, *args):
        """
        Save project_dataset to file directly, unless it hasn't previously been saved.
        :return:
        """
        if not self.project_path:
            self.save_as(None)
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
        if Path.resolve(new_wd) != Path.resolve(self.wd):
            # in fact saving inside a saveas. -> copy from previous dir
            # clear new location, would be overwritten anyway;
            # doing it this roundabout way to keep compatibility with python3.6
            try:
                shutil.rmtree(path=new_wd)
            except FileNotFoundError:
                pass
            shutil.copytree(src=self.wd, dst=new_wd)  # not dirs_exist_ok=True

            # delete old data if it was in the prelim directory
            if self.wd == Path('untitled'):  # Path.cwd() / 'untitled':
                shutil.rmtree(path=self.wd)

            self.wd = new_wd
            # tell the MSA pre-set about it
            if 'aligner' in self.iface:
                gtk_msa.refresh_paths(self)

        self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)
        shared.show_notification(self, 'saved project', stay_secs=2)

    def save_as(self, *args):
        """
        Save with a file dialog. If previously saved, suggest old filename.
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

    def help(self, *args):
        print('HELP!')

    def on_quit(self, *args):
        self.quit()

    def about(self, *args):
        Gtk.AboutDialog(transient_for=self.win, modal=True, title='About AB12PHYLO',
                        authors=['%s<%s>' % (__author__, __email__), 'Dr Remco Stam', 'Corinn Small'],
                        comments='An integrated pipeline for Maximum Likelihood (ML) '
                                 'phylogenetic tree inference from ABI sequencing data',
                        program_name='AB12PHYLO', version=__version__,
                        website='https://gitlab.lrz.de/leokaindl/ab12phylo',
                        website_label='GitLab Repo', license_type=Gtk.License.MIT_X11,
                        logo=GdkPixbuf.Pixbuf.new_from_file(str(ab12phylo_app.ICON))).present()

    def _init_log(self, **kwargs):
        self.log.setLevel(logging.DEBUG)

        if 'filename' in kwargs:
            # init verbose logging to file
            fh = logging.FileHandler(filename=kwargs['filename'], mode='w')
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s\t%(name)s\t%(message)s',
                                              datefmt='%Y-%m-%d %H:%M:%S'))
            self.log.addHandler(fh)

        # init shortened console logging
        sh = logging.StreamHandler(sys.stdout)
        if __verbose__:
            sh.setLevel(logging.DEBUG)
            sh.setFormatter(logging.Formatter('%(levelname)s -- %(message)s'))
        elif __info__:
            sh.setLevel(logging.INFO)
            sh.setFormatter(logging.Formatter('%(message)s'))
        else:
            sh.setLevel(logging.WARNING)

        self.log.addHandler(sh)


def _main():
    app = ab12phylo_app()
    exit_status = app.run(sys.argv)
    sys.exit(exit_status)


if __name__ == '__main__':
    _main()
