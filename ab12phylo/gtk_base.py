# 2021 Leo Kaindl

import ast
import configparser
import hashlib
import logging
import mmap
import pickle
import shutil
import sys
import threading
import warnings
import zipfile
from argparse import Namespace
from pathlib import Path
from time import sleep

import gi
import matplotlib
import pandas as pd
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap

from ab12phylo import repo, ab12phylo_init
from ab12phylo.__init__ import __version__
from ab12phylo.gtk_proj import project_dataset
from ab12phylo_cmd import msa

matplotlib.use('agg')
import matplotlib.pyplot as plt

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, Gio, GLib, GObject, GdkPixbuf

LOG = logging.getLogger(__name__)
__verbose__, __info__ = 0, 1


class ab12phylo_app_base(Gtk.Application):
    TEMPLATE = repo.BASE_DIR / 'ab12phylo' / 'files' / 'gui.glade'
    ICON = repo.BASE_DIR / 'ab12phylo' / 'files' / 'favi.ico'
    CONF = repo.BASE_DIR / 'ab12phylo' / 'conf.cfg'
    CFG = dict()

    def do_activate(self):
        self.add_window(self.win)
        self.win.show_all()

    def do_shutdown(self):
        # delete data if it was in the prelim directory and not saved
        if self.wd == Path('untitled') and not self.project_path:
            LOG.info('shutdown: delete prelim data')
            shutil.rmtree(path=self.wd, ignore_errors=True)

        # kill BLAST thread
        if 'blaster' in self.iface:
            self.iface.pill2kill.set()
            if self.iface.blaster.is_alive():
                self.iface.blaster.stop()
                sleep(1)

        Gtk.Application.do_shutdown(self)

    def do_startup(self):
        Gtk.Application.do_startup(self)
        # connect menu actions and shortcuts
        for act, acc in zip(['new', 'open', '_import', 'save', 'save_as', 'help', 'about',
                             'test', 'on_quit', 'change_colors', 'zoom_in', 'zoom_out'],
                            ['n', 'o', 'i', 's', '<Shift>s', 'h', '<Shift>h',
                             '<Shift>t', 'q', '<Shift>c', 'plus', 'minus']):
            action = Gio.SimpleAction.new(act)
            action.connect('activate', self.__getattribute__(act))  # can pass parameters here
            self.add_action(action)
            self.set_accels_for_action(detailed_action_name='app.%s' % act,
                                       accels=['<Control>%s' % acc])
        self.set_accels_for_action(detailed_action_name='app.help',
                                   accels=['F1', '<Control>h'])

        # fetch the paths from the config
        cfg_parser = configparser.ConfigParser()
        cfg_parser.read(ab12phylo_app_base.CONF)
        if 'Paths' in cfg_parser:
            ab12phylo_app_base.CFG.update(dict(cfg_parser['Paths']))

    def do_command_line(self, cmd):
        options = cmd.get_options_dict().end().unpack()
        if 'version' in options:
            sys.exit('ab12phylo: %s' % __version__)
        if 'initialize' in options or not ab12phylo_app_base.CONF.is_file():
            ab12phylo_init.main()
        if 'open' in options:
            self.load(options['open'])
        if 'proceed' in options:
            self.proceed()
        self.activate()
        return 0

    def __init__(self):
        Gtk.Application.__init__(self, application_id='com.github.lkndl.ab12phylo',
                                 flags=Gio.ApplicationFlags.HANDLES_COMMAND_LINE)
        self.add_main_option('open', ord('o'), GLib.OptionFlags.IN_MAIN,  # ord converts it to integer
                             GLib.OptionArg.STRING, 'To open a project from the commandline',
                             'path to .proj file to open')
        self.add_main_option('proceed', ord('p'), GLib.OptionFlags.IN_MAIN,
                             GLib.OptionArg.NONE, 'Try to proceed to the next stage, do not use')
        self.add_main_option('version', ord('v'), GLib.OptionFlags.IN_MAIN,
                             GLib.OptionArg.NONE, 'Print version information and exit')
        self.add_main_option('initialize', ord('c'), GLib.OptionFlags.IN_MAIN,
                             GLib.OptionArg.NONE, '(re-) download non-python tools and test data. '
                                                  'on Windows, use ab12phylo-init instead.')
        # fetch all named objects from the .glade XML
        iface = dict()
        for widget in Gtk.Builder().new_from_file(str(ab12phylo_app_base.TEMPLATE)).get_objects():
            if widget.find_property('name') and not widget.get_name().startswith('Gtk'):
                iface[widget.get_name()] = widget
        iface = Namespace(**iface)
        self.iface = iface
        self.win = iface.win
        # self.win.set_icon_from_file(str(ab12phylo_app.ICON))

        # set up preliminary working directory
        self.project_path = None
        # self.wd = Path.cwd() / 'untitled'
        self.wd = Path('untitled')  # use relative path should make it movable
        Path.mkdir(self.wd, exist_ok=True)

        self.log = logging.getLogger()
        self._init_log(filename=str(self.wd / 'ab12phylo.log'), mode='w')

        self.win.set_titlebar(iface.tbar)
        self.win.set_hide_titlebar_when_maximized(True)
        self.win.set_title('AB12PHYLO [untitled]')
        LOG.debug('GTK Window initialized')

        # set up an empty thread so Gtk can get used to it
        iface.thread = threading.Thread()
        iface.pill2kill = threading.Event()
        iface.i = 0
        iface.k = 1
        iface.frac = 0
        iface.text = 'idle'
        iface.run_after = list()
        iface.zoomer = Namespace()

        # whether to draw raster (fast) or vector images is a button state

        # get some CSS styling
        mod = b'lighter', b'darker'
        mod2 = 200  # per default, make TreeView text color darker
        theme = Gtk.Settings.get_default().get_property('gtk-theme-name')
        if 'dark' in theme.lower():
            mod = mod[::-1]
            mod2 = 255

        css_provider = Gtk.CssProvider()
        css = b'''
        .seqid { font-size: x-small; padding-left: 3px; padding-right: 2px }
        #gbl_pane separator.wide { min-width: 18px }
        #tree_pane separator { min-width: 3px; }
        progressbar trough progress { 
            min-height: 6px; 
            border-radius: 1px; 
            background-color: @success_color; }
        progressbar trough { min-height: 6px; }
        .hint { opacity: .6; } 
        .prog_labels { font-size: x-small }
        button:active { background-color: #17f }
        notebook header { background-color: transparent; }
        '''
        if 'matcha' in theme.lower():
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

        to_css_color = lambda c: '#' + ''.join([(hex(min(int(c * mod2), 255))[2:]).upper()
                                                for c in list(c)[:-1]])
        iface.FG = to_css_color(sc.get_color(Gtk.StateFlags.ACTIVE))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            iface.BG = to_css_color(sc.get_background_color(Gtk.StateType.NORMAL))

        # prepare shortcuts / accelerators
        self.accelerators = Gtk.AccelGroup()
        self.win.add_accel_group(self.accelerators)

        # connect to the window's delete event to close on x click
        self.win.connect('destroy', lambda *args: self.quit())

        # connect buttons
        iface.next.connect('clicked', self.proceed)
        self.bind_accelerator(self.accelerators, iface.next, '<Alt>Right')
        iface.next.grab_focus()
        iface.back.connect('clicked', self.step_back)
        self.bind_accelerator(self.accelerators, iface.back, '<Alt>Left')
        iface.refresh.connect('clicked', self.re_run)
        # self.bind_accelerator(self.accelerators, iface.refresh, 'Return')
        self.bind_accelerator(self.accelerators, iface.refresh, 'F5')
        iface.reset.connect('clicked', self.reset)
        # any page change
        iface.notebook.connect_after('switch-page', self.refresh)
        # dismissing help
        iface.infobar.connect('response', lambda *args: iface.infobar.set_revealed(False))

        iface.dismiss.connect(
            'clicked', lambda *args: iface.revealer.set_reveal_child(False))
        iface.theme_picker.connect('changed', self.change_theme)
        self.iface.color_dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                            Gtk.STOCK_OK, Gtk.ResponseType.OK)
        self.bind_accelerator(self.accelerators, iface.dismiss, 'Escape')
        self.bind_accelerator(self.accelerators, iface.dismiss, 'Return')

        self.data = project_dataset()
        LOG.debug('interface and dataset initialized')
        # then initialize the notebook pages

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
        Path.mkdir(self.wd, exist_ok=True)
        self.project_path = None
        self.win.set_title('AB12PHYLO [untitled]')
        self.iface.notebook.set_current_page(0)

    def open(self, *args):
        """
        Project file open dialog. Calls self.load internally.
        :return:
        """
        dialog = Gtk.FileChooserDialog(title='open .PROJ project file', parent=None,
                                       action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            self.load(dialog.get_filename())
            dialog.destroy()
            if self.data.page == 7:
                self.refresh()
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
                new_data = pickle.load(proj)  # on win: cannot instantiate 'PosixPath' on your system
            # overwrite content in old dataset in-place rather than re-pointing everything
            self.data.overwrite(new_data)
            try:
                self.init_gene_roll()
            except ValueError as ve:
                LOG.error(ve)
            self.project_path = project_path
            if self.wd == Path('untitled') and project_path.stem != 'untitled':
                shutil.rmtree(self.wd, ignore_errors=True)
            self.wd = self.project_path.parent / self.project_path.stem
            Path.mkdir(self.wd, exist_ok=True)
            self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)
            self.iface.notebook.set_current_page(self.data.page)

            self._init_log(filename=str(self.wd / 'ab12phylo.log'), mode='a')

            if 'df' in self.iface.tempspace:
                del self.iface.tempspace.df
            if self.data.colors:
                repo.colors = list(map(repo.tohex, map(
                    self.data.colors.get, repo.NUCLEOTIDES)))

            # for backwards compatibility: re-name
            if 'raxml_seed' in self.data.ml:
                self.data.ml.ml_seed = self.data.ml.raxml_seed
                del self.data.ml.raxml_seed

        except Exception as e:
            LOG.exception(e)
            self.show_notification('Project could not be loaded', secs=10)

    def save(self, *args, **kwargs):
        """
        Save project_dataset to file directly, unless it hasn't previously been saved.
        :return:
        """
        log_now = LOG.info if 'silent' not in kwargs or not kwargs['silent'] else LOG.debug
        if not self.project_path:
            self.save_as(None)
            return
        with open(self.project_path, 'wb') as proj:
            log_now('saving to %s' % self.project_path)
            try:
                self.data.page = self.iface.notebook.get_current_page()
                pickle.dump(self.data, proj)  # on win: cannot instantiate 'PosixPath' on your system
                log_now('finished save')
            except pickle.PicklingError as pe:
                LOG.error(pe)
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
            try:
                shutil.copytree(src=self.wd, dst=new_wd)  # not dirs_exist_ok=True
            except FileNotFoundError:
                Path.mkdir(new_wd, exist_ok=True)

            if 'blaster' in self.iface and self.iface.blaster.is_alive():
                self.show_message_dialog('Saving to a new location while BLAST is active. '
                                         'Do not delete old data before the search has finished!')
            elif self.wd == Path('untitled'):
                # delete old data if it was in the prelim directory
                try:
                    shutil.rmtree(path=self.wd)
                except FileNotFoundError:
                    pass

            self.wd = new_wd
            log_file = str(self.wd / 'ab12phylo.log')
            try:
                self._init_log(filename=log_file, mode='a')
            except FileNotFoundError:
                LOG.warning(log_file, 'not found')
                self._init_log(filename=log_file, mode='w')

        self.win.set_title('AB12PHYLO [%s]' % self.project_path.stem)
        if 'silent' not in kwargs or not kwargs['silent']:
            self.show_notification('saved project', secs=2)

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
        dialog.set_do_overwrite_confirmation(True)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            self.project_path = Path(dialog.get_filename())
            LOG.debug('got save path to %s' % self.project_path)
            self.save(None)
        dialog.destroy()

    def _import(self, *args):
        """
        Import results from an AB12PHYLO commandline version analysis
        :param args:
        :return:
        """
        dialog = Gtk.FileChooserDialog(title='import commandline project folder',
                                       parent=None, select_multiple=False,
                                       action=Gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        errors = list()
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            folder = Path(dialog.get_filename()).resolve()
            for src, dst in zip(['metadata.tsv', 'msa.fasta', 'tree_TBE.nwk', 'tree_FBP.nwk'],
                                [repo.PATHS.tsv, repo.PATHS.msa, repo.PATHS.tbe, repo.PATHS.fbp]):
                try:
                    shutil.copy(folder / src, self.wd / dst)
                except FileNotFoundError:
                    errors.append('%s not found' % src)
                except Exception as ex:
                    LOG.exception(ex)
            if errors:
                self.show_notification('import failed', errors, secs=10)
                dialog.destroy()
                return

            self.data.genes = False
            self.data.phy.g_lens = False
            for log in ['ab12phylo-p1.log', 'ab12phylo-p2.log', 'ab12phylo-viz.log', 'ab12phylo.log']:
                try:
                    with open(folder / log) as fh:
                        s = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
                        pos = s.find(b'--GENES--')
                        if pos != -1 and not self.data.genes:
                            self.data.genes = s[pos + 10:pos + 190].decode().split('\n')[0].split('::')

                        pos = s.find(b'gene lengths:')
                        if pos != -1 and not self.data.phy.g_lens:
                            self.data.phy.g_lens = {k: v for k, v in ast.literal_eval(
                                s[pos + 14:pos + 190].decode().split('\n')[0])}

                        if self.data.genes and self.data.phy.g_lens:
                            break
                except FileNotFoundError:
                    continue
            if not self.data.genes or not self.data.phy.g_lens:
                self.show_notification('import failed, couldn\'t read '
                                       'annotations from log file', secs=10)
                dialog.destroy()
                return
            LOG.debug('imported from %s' % folder)
            self.iface.notebook.set_current_page(7)
        dialog.destroy()

    def help(self, *args, **kwargs):
        bar = self.iface.infobar
        bar.set_revealed(not bar.get_revealed() or 'show' in kwargs)
        return True

    def test(self, *args):
        if self.wd != Path('untitled'):
            message = 'Open test project, discard unsaved changes to current project?'
            dialog = Gtk.MessageDialog(transient_for=None, flags=0, text=message,
                                       buttons=Gtk.ButtonsType.OK_CANCEL,
                                       message_type=Gtk.MessageType.QUESTION)
            dialog.show_all()  # important
            response = dialog.run()
            dialog.destroy()
            if response != Gtk.ResponseType.OK:
                return

        # make sure to show the infobar
        self.help(show=True)
        # switch directory
        td = Path('test')
        Path.mkdir(td, exist_ok=True)
        with zipfile.ZipFile(repo.BASE_DIR / 'ab12phylo' / 'test_data.zip', 'r') as zf:
            zf.extractall(td / 'test_data')
        shutil.move(td / 'test_data' / 'test.proj', 'test.proj')

        self.load('test.proj')

        # re-build paths
        for model in [self.data.plate_store, self.data.trace_store]:
            for row in model:
                row[0] = str((td / row[0]).resolve())

        # # make relative paths when creating a new test.proj
        # for p in self.data.trace_store:
        #     p[0] = p[0].replace('/home/quirin/PYTHON/AB12PHYLO/ab12phylo/', '')

        # make people not wait:
        self.data.ml.bootstraps = 100
        self.iface.bootstraps.props.text = '100'

    def on_quit(self, *args):
        self.quit()

    def about(self, *args):
        Gtk.AboutDialog(transient_for=self.win, modal=True, title='About AB12PHYLO',
                        authors=['Leo Kaindl<leo.kaindl@tum.de>',
                                 'Remco Stam<stam@wzw.tum.de>', 'Corinn Small'],
                        comments='An integrated pipeline for Maximum Likelihood (ML) '
                                 'phylogenetic tree inference from ABI sequencing data',
                        program_name='AB12PHYLO', version=__version__,
                        website='https://github.com/lkndl/ab12phylo',
                        website_label='GitHub Repo', license_type=Gtk.License.GPL_3_0,
                        logo=GdkPixbuf.Pixbuf.new_from_file(str(ab12phylo_app_base.ICON))).present()

    def _init_log(self, **kwargs):
        self.log.setLevel(logging.DEBUG)

        if 'filename' in kwargs:
            # init verbose logging to file
            fh = logging.FileHandler(filename=kwargs['filename'], mode=kwargs['mode'])
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s\t%(name)s\t%(message)s',
                                              datefmt='%H:%M:%S'))
            self.log.addHandler(fh)

        if 'mode' in kwargs and kwargs['mode'] == 'a':
            return  # otherwise double logging

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

    def refresh(self, *args):
        pass

    def re_run(self, *args):
        pass

    def reset(self, *args):
        pass

    def proceed(self, *args):
        pass

    def step_back(self, *args):
        """Go back to the previous page; set it to un-changed; re-view"""
        self.iface.notebook.prev_page()
        self.data.page = self.iface.notebook.get_current_page()
        self.set_changed(self.data.page, False)
        LOG.debug('stepped back to page %d' % self.data.page)

    # all the shared methods
    def get_errors(self, page):
        return self.data.errors_indicator[page]

    def set_errors(self, page, errors=True):
        self.data.errors_indicator[page] = errors

    def get_changed(self, page):
        return self.data.change_indicator[page]

    def set_changed(self, page, changed=True):
        """
        Set page to changed and disable further stages;
        set page to unchanged and enable next page
        :param gui:
        :param page: the index of the current page
        :param changed: toggle
        :return:
        """
        nb = self.iface.notebook
        indicator = self.data.change_indicator
        if changed:
            # disable later pages
            [page.set_sensitive(False) for page in
             nb.get_children()[page + 1: nb.get_n_pages()]]
            # set later pages to 'changed'
            indicator[page:] = [True] * (nb.get_n_pages() - page)
            indicator[:page] = [False] * page
        else:
            # enable the next page
            nb.get_children()[page + 1].set_sensitive(True)
            # set earlier pages to 'unchanged'
            indicator[:page + 1] = [False] * (page + 1)

    # MARK notify
    @staticmethod
    def show_message_dialog(message, items=None):
        dialog = Gtk.MessageDialog(
            transient_for=None, flags=0, text=message,
            message_type=Gtk.MessageType.WARNING)
        # dialog.format_secondary_text('\n'.join(not_found))
        # not copyable -> not user-friendly
        if items:
            txt_buf = Gtk.TextBuffer()
            txt_buf.props.text = '\n'.join(items)
            txv = Gtk.TextView().new_with_buffer(txt_buf)
            txv.props.wrap_mode = Gtk.WrapMode.WORD_CHAR
            txv.props.margin_end = 20
            txv.props.margin_start = 20
            dialog.get_content_area().add(txv)
        dialog.add_button('OK', 0)
        dialog.get_widget_for_response(0).grab_focus()
        dialog.show_all()  # important
        dialog.run()
        dialog.destroy()

    def show_notification(self, msg, items=None, secs=0):
        revealer = self.iface.revealer
        iface = self.iface
        iface.reveal_title.set_text(msg)
        iface.reveal_list.props.parent.props.visible = items
        if items:
            iface.reveal_list.props.buffer.props.text = '\n'.join(items)
        revealer.set_transition_duration(secs / 4 if secs else 250)
        revealer.set_reveal_child(True)
        if secs:
            threading.Thread(target=self.hide_notification,
                             args=[revealer, secs]).start()

    @staticmethod
    def hide_notification(revealer, secs):
        sleep(secs)
        revealer.set_reveal_child(False)

    # MARK update
    def update(self, page, sensitive=False, pulse=False, *args):
        """
        Keep the progress bar up-to-date, and slowly moving rightwards
        :param page: the index of the current page, which will be frozen
        :param sensitive: allow not freezing the page
        :param pulse: for tasks where no useful duration estimation is possible
        """
        iface = self.iface
        iface.frac = min(
            max(iface.i / iface.k,  # base progress from caller thread iteration
                iface.frac + 0.001),  # pretend to proceed
            (iface.i + 1) / iface.k,  # but do not pass next iteration level
            1)  # and never pass 1
        if iface.thread.is_alive():
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                iface.notebook.get_children()[page].set_sensitive(sensitive)
                if pulse:
                    iface.prog_bar.pulse()
                else:
                    iface.prog_bar.set_fraction(iface.frac)
                for wi in [iface.prog_bar, iface.prog_label]:
                    wi.set_visible(True)
                    wi.set_text(iface.text)
            return True
        else:
            iface.notebook.get_children()[page].set_sensitive(True)
            [wi.set_visible(False) for wi in [iface.prog_bar, iface.prog_label]]
            return False

    def update_ML(self, page, ml, *args):
        """
        Extract the numbers of finished tree searches / bootstrap
        iterations from the RAxML output to update the progress bar.
        """
        iface = self.iface
        if ml.key:
            seen_set = ml.seen[ml.key]
            motif = ml.motifs[ml.key]
            for line in ml.stdout[len(ml.stdout) - 20:]:
                if motif in line:
                    try:
                        seen_set.add(int(line.split(motif)[1].split(',')[0]))
                    except ValueError:
                        pass
            iface.i = len(seen_set) + ml.prev
        return self.update(page, sensitive=True,
                           pulse=(ml.infer_model or ml.ultrafast)
                                 and ml.tool == 'iqtree2')

    def update_BLAST(self, *args):
        """
        Edit the GtkLabel next to the BLAST spinner as a backup progress
        indicator for cases where the spinner is not showing up.
        """
        if 'blast_wrapper' in self.iface and self.iface.blast_wrapper.is_alive():
            spl = self.iface.spin_label
            spl.set_text(' BLAST' + ''.join(['.'] * ((len(spl.get_text()) - 5) % 4)))
            return True
        else:
            self.iface.spin_label.set_text(' BLAST')
            return False

    # MARK TreeViews
    def keep_visible(self, sel, tv, adj, ns):
        """
        For keyboard navigation in previews, scroll the TreeView
        and keep the selection up-to-date
        :param sel: a TreeSelection
        :param adj: the Props of a Vadjustment in a ScrolledWindow
        :param ns: Namespace where previous selection state is stored (amongst others)
        :return:
        """
        mo, tp_iter = sel.get_selected_rows()
        tps = {tp[0] for tp in tp_iter}
        if 'sel' not in ns:
            ns.sel = tps
            return

        try:
            iface = self.iface
            page = iface.notebook.get_current_page()
            h_now = iface.zoomer.sizes[page][-1]
            if h_now * 1.6 < tv.get_allocated_height():
                # do not scroll if it would be jarring
                return
        except KeyError:
            pass

        tp = tps - ns.sel
        ns.sel = tps
        # scrolling only for incremental keyboard selection
        if len(tp) != 1:
            return
        tp = tp.pop()

        if (tp + 1) / len(mo) > (adj.value + adj.page_size) / adj.upper:
            # scroll down
            adj.value = min(adj.upper - adj.page_size,
                            ((tp + 1) / (len(mo) + 1) - adj.page_size / adj.upper) * adj.upper)
        elif tp / len(mo) < adj.value / adj.upper:
            # scroll up
            adj.value = tp / (len(mo) + 1) * adj.upper

    def delete_and_ignore_rows(self, widget, event, page, sel):
        """
        Keep track of the rows that will not be written to the next fasta and delete them from the treeview.
        """
        if Gdk.keyval_name(event.keyval) == 'Delete':
            LOG.debug('delete_and_ignore_rows')
            model, tree_path_iterator = sel.get_selected_rows()
            if page == 2:
                if not self.data.seqdata:
                    self.start_read(run_after=[self.start_trim])
                    return False
                ns = self.data.qal
                if 'ignore_ids' not in ns:
                    ns.ignore_ids = dict()
            elif page == 4:
                ns = self.data.gbl
            else:
                raise RuntimeWarning(f'row deletion called from wrong page {page}')
            for row in reversed(sorted(tree_path_iterator)):
                if page == 2:
                    _id, gene = model[row][:2]
                    _id = self.shift_versions_down_on_deletion(
                        _id, self.data.seqdata[gene], self.data.metadata[gene])
                    if gene in ns.ignore_ids:
                        ns.ignore_ids[gene].add(_id)
                    else:
                        ns.ignore_ids[gene] = {_id}
                    model.remove(model.get_iter(row))
                elif page == 4:
                    ns.ignore_ids.add(model[row][0])
                    model[row][0] = '---'

            self.set_changed(page, True)
            return True
        return False

    @staticmethod
    def shift_versions_down_on_deletion(_id, gene_data, gene_meta):
        match = repo.version_regex.search(_id)
        if match:
            next_v = int(_id[match.start() + 1:]) + 1
            stem = _id[:match.start()]
        else:
            next_v = 1
            stem = _id

        # move entry out of the way
        r = gene_data.pop(_id)
        r.id = _id + '_' + repo.inc_priv_timestamp()
        gene_data[r.id] = r
        gene_meta[r.id] = gene_meta.pop(_id)
        gene_meta[r.id]['quality'] = 'manually dropped at trim_data'

        next_id = '%s.%d' % (stem, next_v)
        while next_id in gene_data:
            # re-organize seqdata and metadata
            r = gene_data.pop(next_id)
            r.id = stem if next_v == 1 else '%s.%d' % (stem, next_v - 1)
            gene_data[r.id] = r
            gene_meta[r.id] = gene_meta.pop(next_id)

            next_v += 1
            next_id = '%s.%d' % (stem, next_v)

        return '%s.%d' % (stem, next_v - 1) if next_v > 1 else _id

    def tv_keypress(self, widget, event, page, selection):
        if Gdk.keyval_name(event.keyval) == 'Delete':
            LOG.debug('delete selected row')
            self.delete_files_from_input_selection(widget, page, selection)
        else:
            LOG.debug('registered keypress, did nothing')

    def delete_files_from_input_selection(self, widget, page, selection, delete_all=False):
        if delete_all:
            delete_all.clear()
            if widget == self.iface.delete_all_trace:
                self.data.trace_store.clear()
        else:
            model, iterator = selection.get_selected_rows()
            [model.remove(model.get_iter(row)) for row in reversed(sorted(iterator))]
        self.set_changed(page, True)
        self.refresh()

    @staticmethod
    def select_seqs(event_box, loc, page, zoom_ns, tv, ns):
        """
        Select sequences from a trim preview directly in the image, in an expected way.
        Shift focus to the labeling treeview on the left, which is already <Delete>-sensitive
        :param event_box:
        :param loc:
        :param sel:
        :param ns:
        :param args:
        :return:
        """
        LOG.debug('select_seqs')
        if tv.props.visible:
            tv.grab_focus()
        sel = tv.get_selection()

        accel_mask = Gtk.accelerator_get_default_mod_mask()
        rect, baseline = event_box.get_allocated_size()
        mo, tree_path_iterator = sel.get_selected_rows()
        idcs = {tp[0] for tp in tree_path_iterator}

        h_now = zoom_ns.sizes[page][3]
        # print('%d:%d:%d:%d' % (tv.get_margin_bottom(), rect.height, h_now, loc.y))
        idx = int((loc.y - (rect.height - h_now) / 2) /
                  h_now / (1 - tv.get_margin_bottom() / tv.get_allocated_height()) * len(mo))
        if idx == '' or idx < 0:
            sel.unselect_all()
            return
        if loc.state & accel_mask == Gdk.ModifierType.CONTROL_MASK:
            sel.select_path(idx) if idx not in idcs else sel.unselect_path(Gtk.TreePath(idx))
        elif loc.state & accel_mask == Gdk.ModifierType.SHIFT_MASK and 'previous' in ns:
            m, n = min(idx, ns.previous), max(idx, ns.previous)
            new_idcs = set(range(m, n))
            if all(idc in idcs for idc in new_idcs):
                sel.unselect_range(Gtk.TreePath(m), Gtk.TreePath(n))
            else:
                sel.select_range(Gtk.TreePath(m), Gtk.TreePath(n))
        else:
            # select only the one clicked on
            sel.unselect_all()
            sel.select_path(idx)
        ns.previous = idx

    @staticmethod
    def save_row_edits(cell, path, new_text, tv, col):
        LOG.debug('save_row_edits')
        mo = tv.get_model()
        old_text = mo[path][col]
        if old_text == new_text:
            return
        mo[path][col] = new_text

    # MARK scaling
    def get_hadj(self):
        return self.iface.zoomer.adj.get_value() * 2

    @staticmethod
    def get_height_resize(widget, event, spacer, scroll_wins, lower=0):
        """
        Adjust the height of a plot in a GtkScrolledWindow depending on the height of the
        associated labeling column. Adjust the width of the spacer below the labels so that
        the scrollbar below the plot has exactly the width of the latter.
        :param widget: the label column whose re-sizing causes this event
        :param event: will be ignored
        :param spacer: below the label column and next to the scrollbar that will be resized, too
        :param scroll_wins: one or two GtkScrolledWindows containing plots
        """
        LOG.debug('re-sizing spacer')
        w, h = widget.get_allocated_width(), widget.get_allocated_height()
        spacer.set_size_request(w, -1)
        if lower:
            h = min(lower, h)
        for sw in scroll_wins:
            LOG.debug('re-sizing %s' % sw.get_name())
            sw.set_max_content_height(h)
            try:
                sw.get_children()[0].set_size_request(w, h)
            except IndexError:
                pass
        return h

    @staticmethod
    def scale(gtk_bin, x, y):
        """
        Stretch a pixbuf residing in a GtkImage, the single child of
        :param gtk_bin:
        :param x: horizontal scaling factor
        :param y: vertical scaling factor
        :return:
        """
        child = gtk_bin.get_child()
        if type(child) != Gtk.Image:
            LOG.info('won\'t rescale %s' % str(type(child)))
            return False
        pb = child.get_pixbuf()
        new_x = min(14000, pb.get_width() * x)
        new_y = min(14000, pb.get_height() * y)
        child.set_from_pixbuf(pb.scale_simple(
            new_x, new_y, GdkPixbuf.InterpType.BILINEAR))
        return new_x, new_y

    def reset_x_scale(self):
        try:
            zoomer = self.iface.zoomer
            with GObject.signal_handler_block(zoomer.adj, zoomer.handle):
                zoomer.adj.props.value = 1
        except Exception as ex:
            LOG.error('re-setting x_scale failed')
            LOG.error(ex)

    def x_scale(self, adj, zoom_ns):
        """Horizontally scale a preview"""
        with GObject.signal_handler_block(adj, zoom_ns.handle):
            a = adj.props
            page = self.iface.notebook.get_current_page()

            a.value = max(.2, a.value)
            x = a.value / zoom_ns.bak
            min_w, min_h, w_now, h_now = zoom_ns.sizes[page]

            if x * w_now > 1.4 * min_w:
                LOG.debug('re-loading images')
                # load larger, so zoom in won't happen so soon again
                a.value = min(a.upper, a.value + 2 * a.step_increment)
                x = a.value / zoom_ns.bak
                if page == 2:
                    self.load_image(zoom_ns, page, self.iface.qal_eventbox,
                                    self.wd / repo.PATHS.preview,
                                    self.data.qal_shape[0] * a.value * 2, h_now)
                elif page == 4:
                    self.load_image(zoom_ns, page, self.iface.gbl_left_vp,
                                    self.wd / repo.PATHS.left,
                                    self.data.msa_shape[0] * a.value * 2, h_now)
                    self.load_image(zoom_ns, page, self.iface.gbl_right_vp,
                                    self.wd / repo.PATHS.right,
                                    self.data.msa_shape[2] * a.value * 2, h_now)
                elif page == 7:
                    self.load_image(zoom_ns, page, self.iface.msa_eventbox,
                                    self.wd / repo.PATHS.phylo_msa,
                                    self.data.phy.shape[0] * a.value * 2, h_now)
                zoom_ns.sizes[page] = [w_now * x, min_h, w_now * x, h_now]
            else:
                LOG.debug('scale x: %.2f fold' % x)
                if page == 2:
                    w_now, h_now = self.scale(self.iface.qal_eventbox, x, 1)
                elif page == 4:
                    self.scale(self.iface.gbl_left_vp, x, 1)
                    w_now, h_now = self.scale(self.iface.gbl_right_vp, x, 1)
                elif page == 7:
                    w_now, h_now = self.scale(self.iface.msa_eventbox, x, 1)
                zoom_ns.sizes[page] = [min(w_now * x, min_w), min_h, w_now * x, h_now]
            zoom_ns.bak = a.value

    def xy_scale(self, widget, event, page):
        """
        Handles zoom in / zoom out on Ctrl+mouse wheel in both x and y
        :param widget:
        :param event:
        :param gui:
        :param page:
        :return:
        """
        iface = self.iface
        accel_mask = Gtk.accelerator_get_default_mod_mask()
        if event.state & accel_mask == Gdk.ModifierType.CONTROL_MASK:
            if not iface.rasterize.props.active:
                LOG.info('won\'t rescale Matplotlib')
                return True
            with GObject.signal_handler_block(iface.zoomer.adj, iface.zoomer.handle):
                self.do_xy_scale(page, event.get_scroll_deltas()[2])
                return True
        return False

    def do_xy_scale(self, page, direction):
        """
        Does the actual zoom in / zoom out on both mouse and key events
        :param page:
        :param direction:
        :return:
        """
        data = self.data
        iface = self.iface
        a = iface.zoomer.adj.props
        bak = a.value
        min_w, min_h, w_now, h_now = iface.zoomer.sizes[page]
        if direction > 0:  # scrolling down -> zoom out -> simple
            # go one tick down
            a.value = max(.2, a.value - a.step_increment)
            if a.value == bak:
                return
            new = a.value / bak
            LOG.debug('scale xy: %.2f fold, %.1f' % (new, a.value))
            if page == 2:
                self.scale(iface.qal_eventbox, new, new)
            elif page == 4:
                self.scale(iface.gbl_left_vp, new, new)
                self.scale(iface.gbl_right_vp, new, new)
            elif page == 7:
                self.scale(iface.msa_eventbox, new, new)
                self.scale(iface.tree_eventbox, new, new)
            # adjust the saved sizes
            iface.zoomer.sizes[page] = [min(min_w * new, min_w), min(min_h * new, min_h), w_now * new,
                                        h_now * new]
        else:
            a.value = min(a.upper, a.value + a.step_increment)
            new = a.value / bak
            if a.value == bak:
                return
            min_w, min_h, w_now, h_now = iface.zoomer.sizes[page]
            if w_now / min_w * new > 3 or h_now / min_h * new > 4:
                # # re-load is due. jump ahead by incrementing again
                # a.value = min(a.upper, a.value + a.step_increment)
                new = a.value / bak
                LOG.debug('re-loading images, scale xy: %.2f fold, %.1f' % (new, a.value))
                if page == 2:
                    h_now = data.qal_shape[1]
                    self.load_image(iface.zoomer, page, iface.qal_eventbox,
                                    self.wd / repo.PATHS.preview,
                                    data.qal_shape[0] * a.value * 2, data.qal_shape[1])
                elif page == 4:
                    h_now = data.gbl_shape[1]
                    self.load_image(iface.zoomer, page, iface.gbl_left_vp,
                                    self.wd / repo.PATHS.left,
                                    data.msa_shape[0] * a.value * 2, data.gbl_shape[1])
                    self.load_image(iface.zoomer, page, iface.gbl_right_vp,
                                    self.wd / repo.PATHS.right,
                                    data.msa_shape[2] * a.value * 2, data.gbl_shape[1])
                elif page == 7:
                    h_now = data.phy.shape[1]
                    self.load_image(iface.zoomer, page, iface.msa_eventbox,
                                    self.wd / repo.PATHS.phylo_msa,
                                    data.phy.shape[0] * a.value * 2, data.phy.shape[1])
                    self.load_image(iface.zoomer, page, iface.tree_eventbox,
                                    self.wd / (data.phy.tx + '.png'), h=data.phy.shape[1])
            else:
                LOG.debug('scale xy: %.2f fold, %.1f' % (new, a.value))
                # scale the easy way
                if page == 2:
                    w_now, h_now = self.scale(iface.qal_eventbox, new, new)
                elif page == 4:
                    self.scale(iface.gbl_left_vp, new, new)
                    w_now, h_now = self.scale(iface.gbl_right_vp, new, new)
                elif page == 7:
                    self.scale(iface.msa_eventbox, new, new)
                    w_now, h_now = self.scale(iface.tree_eventbox, new, new)
            iface.zoomer.sizes[page] = [min_w, min_h, w_now, h_now]
        iface.zoomer.bak = a.value

    def zoom_in(self, *_):
        LOG.debug('zoom in ')
        self.do_xy_scale(self.iface.notebook.get_current_page(), -1)

    def zoom_out(self, *_):
        LOG.debug('zoom out')
        self.do_xy_scale(self.iface.notebook.get_current_page(), 1)

    # MARK images
    @staticmethod
    def load_image(zoom_ns, page, gtk_bin, img_path, w=None, h=None):
        """
        Load an image into a GtkImage child inside gtk_bin, keeping track of observed sizes.
        :param zoom_ns: where all the zooming info is stored, specifically [min_w, min_h,
        w_now, h_now] of the image inside :param gtk_bin inside the zoom_ns.sizes dict.
        """
        # save the current size: only called on enlarging -> easy
        zoom_ns.sizes[page] = [w, h, w, h]

        child = gtk_bin.get_child()
        if child:
            gtk_bin.remove(child)
        child = Gtk.Image()
        gtk_bin.add(child)
        path = str(img_path)
        if h and w:
            pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                filename=path, width=max(1, w), height=h,
                preserve_aspect_ratio=False)
        elif w:
            pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                path, width=w, height=-1, preserve_aspect_ratio=True)
        elif h:
            pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                filename=str(path), width=-1, height=h, preserve_aspect_ratio=True)
        else:
            raise ValueError('Specify at least one out of width / height')
        child.set_from_pixbuf(pb)
        gtk_bin.show_all()

    def load_colorbar(self, gtk_image, gbar=False):
        iface = self.iface
        path = repo.PATHS.gbar if gbar else repo.PATHS.cbar
        try:
            gtk_image.set_from_pixbuf(GdkPixbuf.Pixbuf.new_from_file_at_scale(
                str(self.wd / path), width=250, height=100, preserve_aspect_ratio=True))
        except Exception as ex:
            with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
                LOG.debug('plotting colorbar')
                fig = plt.figure(figsize=(4, .2))
                cax = fig.add_subplot(111)
                i = repo.NUCLEOTIDES.index('-')
                cbar = ColorbarBase(ax=cax, cmap=ListedColormap(repo.colors[:i]),
                                    ticks=[(j + .5) / i for j in range(i)], orientation='horizontal')
                cbar.ax.set_xticklabels(repo.NUCLEOTIDES[:i])
                Path.mkdir(self.wd / path.parent, exist_ok=True)
                fig.savefig(self.wd / path, transparent=True,
                            bbox_inches='tight', pad_inches=0, dpi=600)
                plt.close(fig)
                del fig, cbar
            if gtk_image:
                gtk_image.set_from_pixbuf(GdkPixbuf.Pixbuf.new_from_file_at_scale(
                    str(self.wd / path), width=250, height=100, preserve_aspect_ratio=True))

    # MARK files
    def get_hashes(self, msa_path, page):
        msa_hash = self.file_hash(self.wd / msa_path)
        if msa_hash != self.data.msa_hash:
            self.set_changed(page)
            self.data.msa_hash = msa_hash

    @staticmethod
    def file_hash(file_path):
        h = hashlib.sha256()
        b = bytearray(repo.BUF_SIZE)
        mv = memoryview(b)

        with open(file_path, 'rb', buffering=0) as f:
            for n in iter(lambda: f.readinto(mv), 0):
                h.update(mv[:n])
        return h.hexdigest()

    # MARK numericals
    @staticmethod
    def edit_numerical_entry_up_down(widget, event):
        key = Gdk.keyval_name(event.keyval)
        if key == 'Up':
            widget.set_text(str(1 + int(widget.get_text())))
            return True
        elif key == 'Down':
            widget.set_text(str(max(0, -1 + int(widget.get_text()))))
            return True
        return False

    @staticmethod
    def edit_numerical_entry(widget):
        # filter for numbers only
        value = ''.join([c for c in widget.get_text() if c.isdigit()])
        widget.set_text(value if value else '')  # entering 0 is annoying

    # MARK random
    @staticmethod
    def get_msa_build_cmd(algo, wd, genes, remote=False):
        """
        Initializes a commandline msa_build object
        :param algo:
        :param wd:
        :param genes:
        :param remote:
        :return:
        """
        args = Namespace(**{
            'dir': wd,
            'genes': genes,
            'msa_algo': algo,
            'user': repo.USER,
            'msa': wd / repo.PATHS.msa,
            'sep': repo.SEP,
            'missing_samples': None
        })
        aligner = msa.msa_build(args, None, no_run=True)
        if remote:
            cmd = aligner.build_remote('%s', no_run=True)
        else:
            cmd = aligner.build_local('%s', no_run=True)
        return aligner, cmd

    def write_metadata(self):
        """
        Write the metadata dictionary to the metadata.tsv
        """
        LOG.debug('writing metadata dict to %s' % str(self.wd / repo.PATHS.tsv))
        Path.mkdir(self.wd, exist_ok=True)
        md = self.data.metadata
        # convert metadata dict do pandas DataFrame
        df = pd.concat({gene: pd.DataFrame.from_dict(
            md[gene], orient='index') for gene in md.keys()})
        df.index.names = ['gene', 'id']
        df.reset_index(level=0, inplace=True)
        # write to file
        df.to_csv(self.wd / repo.PATHS.tsv,
                  sep='\t', na_rep='', header=True, index=True)

    @staticmethod
    def bind_accelerator(accelerators, widget, accelerator, signal='clicked'):
        key, mod = Gtk.accelerator_parse(accelerator)
        widget.add_accelerator(signal, accelerators, key, mod, Gtk.AccelFlags.VISIBLE)

    def change_theme(self, combo, *args):
        iface = self.iface
        for wi in iface.flowbox.get_children():
            wi.destroy()
        theme = [repo.technicolor, repo.green_blue, repo.viridis,
                 repo.blue_pink, repo.clustalx][combo.get_active()]
        for nt in repo.NUCLEOTIDES[:6]:
            a = Gtk.ColorButton.new_with_rgba(
                Gdk.RGBA(*theme[nt]))
            a.props.show_editor = True
            iface.flowbox.add(a)
        iface.flowbox.show_all()
        LOG.debug('changed color theme')

    def change_colors(self, *args):
        iface = self.iface
        if not self.data.colors:
            self.data.colors = {nt: c for nt, c in repo.technicolor.items()}

        for wi in iface.flowbox.get_children():
            wi.destroy()

        for nt in repo.NUCLEOTIDES[:6]:
            a = Gtk.ColorButton.new_with_rgba(
                Gdk.RGBA(*self.data.colors[nt]))
            a.props.show_editor = True
            iface.flowbox.add(a)
        iface.flowbox.show_all()

        response = self.iface.color_dialog.run()
        if response == Gtk.ResponseType.OK:
            for wi, nt in zip(iface.flowbox.get_children(), repo.NUCLEOTIDES[:6]):
                self.data.colors[nt] = tuple(wi.get_child().get_rgba())

            repo.colors = list(map(repo.tohex, map(self.data.colors.get, repo.NUCLEOTIDES)))
            LOG.debug('changed colors')
            self.load_colorbar(None)  # force a re-plot
        self.iface.color_dialog.hide()
        self.re_run()

