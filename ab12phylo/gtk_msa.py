# 2021 Leo Kaindl

import logging
import shutil
import subprocess
import sys
import threading
from pathlib import Path

import gi
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import GObject

from ab12phylo import repo
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 3


class msa_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.msa_algo.set_entry_text_column(0)
        iface.msa_algo.set_id_column(0)
        iface.remote_algo.set_id_column(0)

        iface.msa_cmd.connect('focus_out_event', lambda widget, *args: data.msa.cmd.update(
            {repo.toalgo(iface.msa_algo.get_active_text()): widget.get_buffer().props.text.strip()}))
        iface.remote_cmd.connect('focus_out_event', lambda widget, *args: data.msa.remote_cmd.update(
            {repo.toalgo(iface.remote_algo.get_active_text()): widget.get_buffer().props.text.strip()}))

        iface.msa_algo.connect('changed', self.get_help)
        iface.remote_algo.connect('changed', self.get_help, True)
        iface.msa_import.connect('file-set', self.load_msa)
        iface.msa_exe.connect('file-set', self.get_help, True, True)

        # connect buttons
        iface.msa_build.connect('clicked', self.start_align)
        self.bind_accelerator(self.accelerators, iface.msa_build, 'Return')
        iface.remote_build.connect('clicked', self.start_align, True)
        self.bind_accelerator(self.accelerators, iface.remote_build, 'Return')

        data.msa.stack_child_name = iface.align_stack.get_visible_child_name()

    def get_help(self, widget, remote=False, try_path=False):
        data = self.data
        iface = self.iface

        if not data.genes:
            return
        if remote:
            data.msa.algo = repo.toalgo(iface.remote_algo.get_active_text())
            client = repo.TOOLS / 'MSA_clients' / (data.msa.algo + '.py')
            # if whitespace in the python path, encapsulate with ""
            py3 = sys.executable
            py3 = py3 if ' ' not in py3 else f'"{py3}"'

            self.set_helpers(f'{py3} {client} ', iface.remote_help,
                             data.msa.remote_cmd, data.msa.algo, True, iface.remote_cmd)
        else:
            data.msa.algo = repo.toalgo(iface.msa_algo.get_active_text())
            exe = shutil.which(data.msa.algo)
            exe = widget.get_active_text() if try_path else exe
            if exe:
                # get the --help output and save it in the lookup field on the right
                iface.msa_exe.set_filename(exe)
                self.set_helpers('%s --help; exit 0' % exe, iface.msa_help,
                                 data.msa.cmd, data.msa.algo, False, iface.msa_cmd)
            else:
                # no executable found; unselect in path box
                iface.msa_exe.unselect_all()
                txt = data.msa.algo + ' was not found on your system $PATH. You can try ' \
                                      'manually specifying the path to the executable.'
                iface.msa_help.get_buffer().props.text = txt  # show this snarky line above
                iface.msa_cmd.get_buffer().props.text = ''  # no cmd suggestion

    def set_helpers(self, cmdline, help_view, help_dict, algo, remote, cmd_view):
        proc = subprocess.Popen(cmdline, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        txt = out.decode().strip() + '\n\n' + err.decode().strip()
        help_view.get_buffer().props.text = txt.strip()

        # get the suggested command and allow user modification in the left field
        txt = help_dict.get(algo, '')  # fetch saved
        if txt == '' or 'aligner' not in self.iface:  # deleting all content will also get you back the original
            self.iface.aligner, txt = self.get_msa_build_cmd(algo, self.wd, self.data.genes, remote)  # fetch new
            help_dict[algo] = txt  # save

        cmd_view.get_buffer().props.text = help_dict[algo]  # show

    def refresh_paths(self):
        LOG.debug('refresh_paths')
        old = str(self.iface.aligner.dir)
        self.iface.aligner.reset_paths(self.wd, self.wd / repo.PATHS.msa)
        for algo, txt in self.data.msa.cmd.items():
            self.data.msa.cmd[algo] = txt.replace(old, str(self.wd))
        # also replace it in the GtkTextView
        bf = self.iface.msa_cmd.get_buffer()
        bf.props.text = bf.props.text.replace(old, str(self.wd))

    def start_align(self, widget, remote=False, run_after=None):
        """
        Starts an MSA building thread unless one of the following conditions is met:
        a) another thread is active -> abort + forbid proceeding.
        b) this function was called from the _Next button and an MSA already exists
        -> accept + allow proceeding.
        c) there were no changes registered for this page and no proceeding (not
        called by _Next) -> accept + show notification.

        :param widget: required for callback, ignored
        :param remote: if the MSA should be constructed using the EMBL-EBI API
        :param run_after: the function to run afterwards; usually flip to next page
        :return:
        """
        data = self.data
        iface = self.iface
        if iface.thread.is_alive():  # a)
            self.show_notification('Busy', secs=1)
            return
        elif data.msa.last_cmd == [data.msa.cmd, data.msa.remote_cmd][remote][data.msa.algo] \
                and not self.get_changed(PAGE):
            self.show_notification('MSA already generated, please proceed')
            return
        elif all([(self.wd / gene / ('%s_raw_msa.fasta' % gene)).exists() for gene in data.genes]) \
                and run_after and not self.get_errors(PAGE):  # b)
            self.set_changed(PAGE, False)
            [do_func() for do_func in run_after]
            return
        if 'aligner' not in iface:
            self.get_help(None, remote)
        self.save_ui_state()
        self.save()

        data.msa_lens.clear()
        iface.align_stack.props.sensitive = False
        iface.thread = threading.Thread(target=self.do_align, args=[remote])
        iface.run_after = run_after
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        return
        # return to main loop

    def do_align(self, remote=False):
        data = self.data
        iface = self.iface
        errors = list()
        iface.frac = .05
        iface.i = 0
        iface.k = len(data.genes)
        funcs, arg_dicts = [iface.aligner.build_local, iface.aligner.build_remote], \
                           [data.msa.cmd, data.msa.remote_cmd]
        try:
            for gene in data.genes:
                iface.text = 'aligning %s [%d/%d]' % (gene, iface.i + 1, iface.k)
                LOG.debug(iface.text)
                arg = arg_dicts[remote][data.msa.algo]
                data.msa.last_cmd = arg
                try:
                    funcs[remote](gene, new_arg=
                    arg % tuple([gene] * (4 - remote)))  # interpreting bool as int here
                except (FileNotFoundError, subprocess.CalledProcessError):
                    self.refresh_paths()
                    # try again once more
                    funcs[remote](gene, new_arg=arg % tuple([gene] * 4))
                # fetch MSA length
                for r in SeqIO.parse(self.wd / gene / ('%s_raw_msa.fasta' % gene), 'fasta'):
                    data.msa_lens.append(len(r))
                    break
                iface.i += 1
            iface.frac = 1
            iface.text = 'idle'
        except (OSError, subprocess.CalledProcessError) as e:
            errors.append('%s at task %d (%s). invalid command?'
                          % (type(e), iface.i, iface.text))
        except FileNotFoundError:
            errors.append('MSA/sequences file not found. Did you just save somewhere new?')
        except TypeError:
            raise ValueError('replace string formatting with something smarter')
        GObject.idle_add(self.stop_align, errors)

    def stop_align(self, errors):
        iface = self.iface
        iface.thread.join()
        iface.align_stack.props.sensitive = True
        self.win.show_all()
        LOG.info('msa thread idle')
        self.set_errors(PAGE, bool(errors))
        self.set_changed(PAGE, False)
        if errors:
            self.show_notification('Errors during MSA building', errors)
            return
        if iface.run_after:
            [do_func() for do_func in iface.run_after]
        else:
            self.show_notification('MSA building finished', items=None)
        return

    def load_msa(self, widget):
        data = self.data
        iface = self.iface
        try:
            Path.mkdir(self.wd / repo.PATHS.import_msa.parent, exist_ok=True)
            shutil.copy(widget.get_filename(), self.wd / repo.PATHS.import_msa)
        except shutil.SameFileError:
            pass
        except Exception as ex:
            self.show_notification(str(ex))
            LOG.error(ex)
        self.get_hashes(repo.PATHS.import_msa, PAGE)
        data.genes = ['import']
        data.gene_ids = {'import': {r.id for r in SeqIO.parse(
            self.wd / repo.PATHS.import_msa, 'fasta')}}
        iface.aligner, cmd = self.get_msa_build_cmd(
            repo.toalgo(iface.msa_algo.get_active_text()), self.wd, data.genes)
        # write a metadata.tsv
        with open(self.wd / repo.PATHS.tsv, 'w') as metadata:
            metadata.write('id\tgene\n')
            for _id in data.gene_ids['import']:
                metadata.write('%s\timport\n' % _id)
        LOG.debug('using imported MSA')
        self.set_changed(PAGE, False)
        self.save_ui_state()

    def refresh(self):
        if self.iface.align_stack.get_visible_child_name() == 'local':
            self.get_help(None)
        elif self.iface.align_stack.get_visible_child_name() == 'remote':
            self.get_help(None, remote=True)

    def save_ui_state(self):
        ns = self.data.msa
        iface = self.iface
        ns.stack_child_name = iface.align_stack.get_visible_child_name()
        ns.msa_algo_id = iface.msa_algo.get_active_id()
        ns.msa_exe_filename = iface.msa_exe.get_filename()
        ns.remote_algo_id = iface.remote_algo.get_active_id()
        ns.msa_import_filename = iface.msa_import.get_filename()

    def reload_ui_state(self):
        ns = self.data.msa
        iface = self.iface
        iface.align_stack.set_visible_child_name(ns.stack_child_name)
        iface.msa_algo.set_active_id(ns.msa_algo_id)
        if ns.msa_exe_filename:
            iface.msa_exe.set_filename(ns.msa_exe_filename)
        iface.remote_algo.set_active_id(ns.remote_algo_id)
        if ns.msa_import_filename:
            iface.msa_import.set_filename(ns.msa_import_filename)

        iface.msa_cmd.get_buffer().props.text = \
            ns.cmd.get(repo.toalgo(ns.msa_algo_id), '')
        iface.remote_cmd.get_buffer().props.text = \
            ns.remote_cmd.get(repo.toalgo(ns.remote_algo_id), '')
