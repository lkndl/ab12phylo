# 2021 Leo Kaindl

import logging
import mmap
import random
import shlex
import shutil
import stat
import sys
import threading
from os import cpu_count
from pathlib import Path
from subprocess import call, run, Popen, PIPE
from time import sleep, time
from zipfile import ZipFile, ZIP_DEFLATED

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo import repo
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 6
CREATE_NO_WINDOW = 0x8000000


class ml_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.ml_tool.set_entry_text_column(0)
        iface.ml_tool.set_id_column(0)
        iface.ml_tool.connect('changed', self._load_tool_help)
        iface.ml_tool.connect_after('changed', lambda widget, *args: self.start_ML(widget, 'prep'))
        iface.ml_exe.connect('file-set', self._load_tool_help, True)
        iface.ultrafast.connect('toggled', self._ultrafast_hint)
        for wi in [iface.ultrafast, iface.infer_model]:
            wi.connect_after('toggled', lambda widget, *args: self.start_ML(widget, 'prep'))

        # connect the parser to basically everything
        for wi in [iface.bootstraps, iface.ml_seed, iface.cpu_use, iface.rand, iface.pars,
                   iface.evo_modify, iface.evo_model, iface.ml_help, iface.ml_cmd]:
            wi.connect('focus_out_event', self._parse_ml)

        iface.evo_model.set_model(data.evo_models)
        iface.evo_model.get_child().connect(
            'focus_out_event', self._change_evo_model, iface.evo_modify, data.ml)
        iface.evo_block = iface.evo_model.connect_after(
            'changed', self._load_model_file, iface.evo_modify)
        iface.evo_model.handler_block(iface.evo_block)

        iface.ml_run.connect('clicked', self.start_ML, 'run')
        iface.ml_export.connect('clicked', self.start_ML, 'export')
        iface.ml_import.connect('clicked', self.import_tree)

        for wi in [iface.bootstraps, iface.rand, iface.pars, iface.ml_seed]:
            wi.connect('key-press-event', self.edit_numerical_entry_up_down)
            wi.connect('changed', self.edit_numerical_entry)

        iface.ml_page_seen = False

    def _ultrafast_hint(self, widget, *args):
        """
        If the ultrafast bootstrapping button is toggled on and the current
        number of bootstraps is lower than 1000, hint that this won't work.
        """
        if not widget.get_active():
            return
        wi = self.iface.bootstraps
        val = int([i for i in [wi.get_text(), wi.get_placeholder_text()] if i][0])
        if val >= 1000:
            return
        self.show_notification(msg='ultrafast bootstrapping will only work '
                                   'with ≥ 1000 iterations', secs=10)

    def _parse_ml(self, widget, event_focus):
        """An interface parser to save changes to the UI in the model."""
        data = self.data
        iface = self.iface
        ml = data.ml
        ml.binary = iface.ml_exe.get_filename()
        for w_name in ['evo_modify', 'bootstraps', 'rand', 'pars', 'ml_seed', 'cpu_use']:
            wi = iface.__getattribute__(w_name)
            val = [i for i in [wi.get_text(), wi.get_placeholder_text()] if i][0]
            if w_name in ['bootstraps', 'rand', 'pars']:
                val = int(val)
            ml.__setattr__(w_name, val)
        ml.evo_model = data.evo_models[iface.evo_model.get_active()]
        ml.infer_model = iface.infer_model.get_active()
        ml.ultrafast = iface.ultrafast.get_active()
        tx = iface.ml_cmd.get_buffer().props.text
        if tx:
            ml.ml_cmd = tx
        if ml.evo_model[1]:
            # do not allow modifier suffixes if a model file is used
            ml.evo_modify = ''
            ml.evo_model = str(self.wd / 'RAxML' / 'user_model')
        else:
            ml.evo_model = ml.evo_model[0]

    def import_tree(self, widget):
        """Import a tree file or two"""
        dialog = Gtk.FileChooserDialog(title='import tree(s)',
                                       parent=None, select_multiple=True,
                                       action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            paths = [Path(p).resolve() for p in dialog.get_filenames()]
            if len(paths) > 2:
                self.show_message_dialog('Please select at most two tree files.')
            else:
                errors = list()
                for path in paths:
                    try:
                        if 'FBP' in path.name.upper():
                            shutil.copy(path, self.wd / repo.PATHS.fbp)
                        elif 'TBE' in path.name.upper():
                            shutil.copy(path, self.wd / repo.PATHS.tbe)
                        else:
                            errors.append(path.name)
                    except shutil.SameFileError as sfe:
                        errors.append(f'{path} is already in the project')
                if errors:
                    self.show_message_dialog('Not immediately recognized as either '
                                             'tree_FBP.nwk or tree_TBE.nwk. You can also copy '
                                             'it/them to %s manually.' % self.wd, items=errors)
                else:
                    self.show_notification('imported trees:',
                                           [path.name for path in paths], 2)
                    self.set_errors(PAGE, False)

        dialog.destroy()
        self.set_changed(PAGE, False)

    @staticmethod
    def _change_evo_model(entry, event_focus, evo_modify, ml):
        combo = entry.get_parent().get_parent()
        if combo.get_active_iter():
            ml.evo_model = entry.get_text()
            evo_modify.props.sensitive = True
        else:
            tx = entry.get_text()
            combo.get_model().append([tx, None])
            ml.evo_model = tx  # save in project dataset
            LOG.debug('entered custom evo model %s' % tx)

    def _load_model_file(self, combo, evo_modify):
        if combo.get_active_iter() and combo.get_active_id() == 'from file':
            dialog = Gtk.FileChooserDialog(title='select file with partition table',
                                           parent=None, select_multiple=False,
                                           action=Gtk.FileChooserAction.OPEN)
            dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                               Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = dialog.run()
            if response == Gtk.ResponseType.OK:
                tx = Path(dialog.get_filename()).resolve()
                try:
                    Path.mkdir(self.wd / 'RAxML', exist_ok=True)
                    shutil.copy(tx, self.wd / 'RAxML' / 'user_model')
                    shutil.copy(tx, self.wd / 'RAxML' / 'ml.raxml.bestModel')
                except FileNotFoundError:
                    self.show_notification('File does not exist', secs=1)
                    dialog.destroy()
                    return
                except shutil.SameFileError:
                    pass
                evo_modify.props.sensitive = False
                combo.get_model().append([tx.name, str(tx)])
                combo.set_active_id(tx.name)
                LOG.debug('selected partitioned model file: %s -> RAxML/user_model' % str(tx))
            dialog.destroy()

    def _load_tool_help(self, widget, try_path=False):
        """
        Load the RAxML-NG / IQ-Tree path from the config file or use the user-defined one.
        If --help exits with OSError, despair.
        """
        iface = self.iface
        ml = self.data.ml
        LOG.debug(f'_load_tool_help {repo.toname(ml.tool)}')

        if not iface.ml_page_seen:
            iface.ml_tool.set_active_id(repo.toname(ml.tool))
        else:
            ml.tool = repo.toalgo(iface.ml_tool.get_active_text())

        for wi in [iface.iqtree_label, iface.ultrafast, iface.infer_model]:
            wi.set_sensitive(ml.tool != 'raxml-ng')

        if try_path:
            binary = widget.get_filename()
        else:
            binary = ab12phylo_app_base.CFG.get(ml.tool, shutil.which(ml.tool))

        if binary:
            # May be missing, for example no RAxML-NG on Windows
            ml.binary = binary
            iface.ml_exe.set_filename(binary)
            binary = Path(binary)
            # Ensure the file is executable, and present
            try:
                binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
            except FileNotFoundError:
                LOG.error(f"File for {iface.ml_tool.get_active_text()} not found "
                          f"at {binary}. Please re-run 'ab12phylo --initialize' or "
                          f"manually correct the path in {ab12phylo_app_base.CONF}")

            try:
                res = run(args=f'{ml.binary} --help', shell=True,
                          stdout=PIPE, stderr=PIPE)
                iface.ml_help.props.buffer.props.text = res.stdout.decode().lstrip()
                LOG.debug('got %s --help' % ml.tool)
            except OSError as osex:
                LOG.exception(osex)
                if sys.platform in ['win32', 'darwin', 'cygwin']:
                    self.show_notification('The selected binary doesn\'t work for this OS.\n'
                                           'Please use a different tool or export a .zip\n'
                                           'and run ML inference on a different machine.\n'
                                           + ml.binary)
                else:
                    assert sys.platform == 'linux', 'What\'s AIX?'

        try:
            # get number of available CPUs
            cpus = cpu_count()
            iface.cpu_count.set_text(str(cpus))
            cpu_adj = iface.cpu_use.get_adjustment().props
            cpu_adj.upper = cpus
            cpu_adj.value = cpus
            LOG.debug('found %d CPUs' % cpus)
        except Exception as ex:
            LOG.error('reading CPUs failed')

    def reload_ui_state(self):
        data = self.data
        iface = self.iface
        ml = data.ml
        for w_name in ['bootstraps', 'rand', 'pars']:
            iface.__getattribute__(w_name).set_text(str(ml.__getattribute__(w_name)))
        iface.ml_seed.set_text(str(ml.ml_seed) if 'ml_seed' in ml else '')
        iface.evo_model.set_active_id(ml.evo_model)
        iface.evo_modify.set_text(ml.evo_modify)
        iface.in_shell.set_active(ml.in_shell)
        iface.infer_model.set_active(ml.infer_model)
        iface.ultrafast.set_active(ml.ultrafast)

        iface.ml_tool.set_active_id(repo.toname(ml.tool))
        iface.ml_cmd.get_buffer().props.text = ml.ml_cmd
        for wi in [iface.iqtree_label, iface.ultrafast, iface.infer_model]:
            wi.set_sensitive(ml.tool != 'raxml-ng')

    def refresh(self):
        """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
        LOG.debug('ML refresh')
        data = self.data
        iface = self.iface
        if not iface.ml_page_seen:
            self._load_tool_help(None)
            self.reload_ui_state(PAGE)
            iface.evo_model.handler_unblock(iface.evo_block)
            iface.evo_model.set_active_id(data.ml.evo_model)
            if sys.platform in ['win32', 'cygwin']:
                iface.in_shell.set_active(False)
                iface.in_shell.set_sensitive(False)
            iface.ml_page_seen = True
            self.start_ML(None, 'prep')

        self.set_errors(PAGE, not any((self.wd / a).is_file()
                                      for a in [repo.PATHS.tbe, repo.PATHS.fbp]))
        # change the button to its default state
        (im, la), tx = iface.ml_run.get_child().get_children(), 'Run'
        im.set_from_icon_name('media-playback-start-symbolic', 4)
        la.set_text('Run')

    def start_ML(self, widget, mode, run_after=None):
        """
        Set-up the ML inference thread.
        :param widget: ignored, to enable use as callback
        :param mode: from {'prep', 'run', 'export'}
        :param run_after: for chaining calls, run after thread finishes
        """
        data = self.data
        iface = self.iface
        ml = data.ml

        if not data.genes:
            self.show_notification('No genes', secs=1)
            return
        elif iface.thread.is_alive() and mode == 'prep':
            return
        elif iface.thread.is_alive():
            # The button was pressed when it was in the 'Stop' state
            iface.pill2kill.set()
            return

        # parse interface
        # ml.tool = repo.toalgo(iface.ml_tool.get_active_text())
        self._parse_ml(None, None)

        ml.ml_seed = random.randint(0, max(1000, ml.bootstraps)) \
            if ml.ml_seed == 'random' else int(ml.ml_seed)
        iface.ml_seed.props.text = str(ml.ml_seed)
        Path.mkdir(self.wd / 'RAxML', exist_ok=True)
        Path.mkdir(self.wd / 'IQ-Tree', exist_ok=True)
        ml.in_shell = iface.in_shell.get_active()

        iface.run_after = run_after
        if ml.infer_model:
            iface.run_after = [self._unset_modelfinder]
        # to keep the progress bar up-to-date:
        ml.prev = 0
        ml.key = False
        ml.stdout = list()
        ml.seen = {'ML': set(), 'BS': set()}
        ml.motifs = {'ML': 'ML tree search #', 'BS': 'Bootstrap tree #'}
        iface.k = ml.bootstraps + ml.rand + ml.pars + 3 if mode == 'run' else 2
        if mode != 'prep':
            self.save()
            sleep(.1)

            # change the button to its stop state
            (im, la), tx = iface.ml_run.get_child().get_children(), 'Run'
            im.set_from_icon_name('media-playback-stop-symbolic', 4)
            la.set_text('Stop')
            iface.tup = (im, la, tx)

        iface.thread = threading.Thread(target=self.do_ML, args=[mode])
        GObject.timeout_add(100, self.update_ML, PAGE, ml)
        iface.thread.start()
        return  # to main loop

    def do_ML(self, mode):
        """Run the ML inference thread"""
        data = self.data
        iface = self.iface
        ml = data.ml
        start = time()
        msa = self.wd / repo.PATHS.msa
        prefix = self.wd / ('RAxML' if ml.tool == 'raxml-ng' else 'IQ-Tree')
        shell = self.wd / f"{'raxml' if ml.tool == 'raxml-ng' else 'iqtree'}_run.sh"
        errors = list()
        msg = False
        iface.i = 0
        assert ml.tool in {'raxml-ng', 'iqtree2'}, f'got unexpected tool: {ml.tool}'

        if mode == 'prep' or (mode == 'export' and ml.infer_model and ml.tool == 'iqtree2'):
            # Create the calls and write to the TextBuffer
            if ml.tool == 'raxml-ng':
                chck = '# check MSA\n' \
                       '"{raxml-ng}" --msa "{msa}" --check --model ' \
                       '"{evo_model}{model_modifier}" --prefix "{chk_prefix}"'

                inML = '# infer ML tree\n' \
                       '"{raxml-ng}" --msa "{msa}" --model ' \
                       '"{evo_model}{model_modifier}" --prefix "{ml_prefix}" ' \
                       '--seed {seed} --threads auto{{{cpus}}} --workers auto{{{cpus}}} ' \
                       '--redo --tree {start_trees}'

                boot = '# bootstrapping\n' \
                       '"{raxml-ng}" --bootstrap --msa "{msa}" --model "{inferred_model}" ' \
                       '--tree "{ml_tree}" --prefix "{bs_prefix}" --bs-trees {bootstraps} ' \
                       '--seed {seed} --threads auto{{{cpus}}} --workers auto{{{cpus}}} --redo'

                supp = '# calc. branch support\n' \
                       '"{raxml-ng}" --support --tree "{ml_tree}" --bs-trees "{bs_trees_file}" ' \
                       '--bs-metric fbp,tbe --prefix "{sp_prefix}" ' \
                       '--threads auto{{{cpus}}} --workers auto{{{cpus}}} --redo'

                iface.ml_cmd.get_buffer().props.text = '\n\n'.join((chck, inML, boot, supp))
            elif ml.tool == 'iqtree2':
                calls = list()

                if ml.infer_model and mode != 'export':
                    calls.append('# find best model\n'
                                 '"{iqtree2}" -s "{msa}" -m MF -mtree -mset raxml '
                                 '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                                 '--prefix "{mf_prefix}" ')
                else:
                    calls.append('# infer ML tree\n'
                                 '"{iqtree2}" -s "{msa}" -m "{evo_model}{model_modifier}" '
                                 '-ninit {start_trees} -ntop {start_trees} '
                                 '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                                 '--prefix "{ml_prefix}" ')
                    if not ml.ultrafast:
                        calls.append('# non-parametric bootstrapping\n'
                                     '"{iqtree2}" -s "{msa}" -m "{evo_model}{model_modifier}" '
                                     '-te "{ml_prefix}.treefile" -b {bootstraps} '
                                     '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                                     '--prefix "{bs_prefix}" -quiet')
                    else:
                        calls.append('# ultrafast bootstrapping\n'
                                     '"{iqtree2}" -s "{msa}" -m "{evo_model}{model_modifier}" '
                                     '-t "{ml_prefix}.treefile" -B {bootstraps} -wbtl '
                                     '--seed {seed} -nt AUTO -ntmax {cpus} -redo '
                                     '--prefix "{uf_prefix}" ')
                    calls.append('# calc. TBE branch support\n'
                                 '"{iqtree2}" -sup "{ml_prefix}.treefile" '
                                 '-t "{boot_trees}" --tbe '
                                 '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                                 '--prefix "{sp_prefix}" ')

                iface.ml_cmd.get_buffer().props.text = '\n\n'.join(calls)

            if mode == 'prep':
                GObject.idle_add(self.stop_prep, errors)
                return True

        assert mode in {'run', 'export'}, f'got unexpected mode: {mode}'

        # Read call templates from buffer
        calls, descs, limits, keys, format_args = list(), list(), list(), list(), dict()
        ml.ml_cmd = iface.ml_cmd.get_buffer().props.text
        for line in ml.ml_cmd.strip().split('\n'):
            if line.startswith('#'):
                descs.append(line.strip()[2:])
            elif line:
                calls.append(line)

        # Prepare args to fill call templates
        if ml.tool == 'raxml-ng':
            keys = [False, 'ML', 'BS', False]
            limits = [0, 1, ml.rand + ml.pars + 1, iface.k - 2]
            format_args = {'raxml-ng': ml.binary,
                           'msa': msa, 'seed': ml.ml_seed,
                           'evo_model': ml.evo_model,
                           'model_modifier': ml.evo_modify,
                           'chk_prefix': prefix / 'chk',
                           'ml_prefix': prefix / 'ml',
                           'bs_prefix': prefix / 'bs',
                           'sp_prefix': prefix / 'sp',
                           'cpus': ml.cpu_use,
                           'bootstraps': ml.bootstraps,
                           'inferred_model': prefix / 'ml.raxml.bestModel',
                           'ml_tree': prefix / 'ml.raxml.bestTree',
                           'bs_trees_file': prefix / 'bs.raxml.bootstraps',
                           'start_trees': ','.join(
                               [a for a in ['rand{%d}' % ml.rand if ml.rand > 0 else None,
                                            'pars{%d}' % ml.pars if ml.pars > 0 else None] if a])}
        elif ml.tool == 'iqtree2':
            keys = ['ML', 'BS', False]
            limits = [1, ml.rand + ml.pars, iface.k - 2]
            format_args = {'iqtree2': ml.binary,
                           'msa': msa, 'seed': ml.ml_seed,
                           'evo_model': ml.evo_model,
                           'model_modifier': ml.evo_modify,
                           'mf_prefix': prefix / 'mf',
                           'ml_prefix': prefix / 'ml',
                           'bs_prefix': prefix / 'bs',
                           'uf_prefix': prefix / 'uf',
                           'sp_prefix': prefix / 'sp',
                           'cpus': ml.cpu_use,
                           'bootstraps': ml.bootstraps,
                           'boot_trees': prefix / ('uf.ufboot' if ml.ultrafast else 'bs.boottrees'),
                           'start_trees': ml.rand + ml.pars}

        if ml.in_shell:
            with open(shell, 'w') as sh:
                sh.write('#!/bin/bash\n\n')

        # loop over the stages
        for i, (arg, desc, key, prev) in enumerate(zip(calls, descs, keys, limits)):

            if ml.in_shell and mode == 'run':
                if ml.tool == 'raxml-ng':
                    with open(shell, 'a') as sh:
                        sh.write(f'# {desc}\n')
                        # sh.write(arg.format(**format_args)
                        if i != 2:
                            sh.write(arg.format(**format_args))
                        else:
                            # special case bootstrapping:
                            fifo = Path('pipe')
                            if fifo.exists():
                                fifo.unlink()
                            bash = '''
mkfifo pipe || exit 1
(%s) > pipe &
pid=$!
echo "AB12PHYLO: Bootstrapping PID is $pid"
while read -r line; do
    echo "$line"
    if [[ "${line::12}" == "Elapsed time" ]]; then
        echo "AB12PHYLO: Finished bootstrapping, terminating process $pid to ensure it exits."
        kill -s SIGTERM $pid
        break
    fi
done < pipe
rm pipe
                            '''.strip() % (arg.format(**format_args))
                            sh.write(bash)
                        sh.write('\n\necho "AB12PHYLO: %s done"\n' % desc)
                        if i != 3:
                            sh.write('\nsleep 1s\n\n')
                        continue

                elif ml.tool == 'iqtree2':
                    with open(shell, 'a') as sh:
                        sh.write(f'# {desc}\n')
                        sh.write(arg.format(**format_args))
                        sh.write('\n\necho "AB12PHYLO: %s done"\n' % desc)
                        continue

            # running RAxML or IQ-Tree, but live rather than in shell mode
            iface.text = desc
            LOG.info(iface.text)
            ml.stdout = list()
            ml.key = key
            ml.prev = prev

            # For RAxML, the first call is the MSA check - which is quick and useful.
            if not (mode == 'export' and ml.tool == 'iqtree2'):
                # read realtime RAxML output line by line
                if sys.platform in ['win32', 'cygwin']:
                    proc = Popen(args=shlex.split(arg.format(**format_args)),
                                 stdout=PIPE, stderr=PIPE, creationflags=CREATE_NO_WINDOW)
                else:
                    proc = Popen(args=shlex.split(arg.format(**format_args)),
                                 stdout=PIPE, stderr=PIPE)
                while True and not iface.pill2kill.is_set():
                    line = proc.stdout.readline()
                    if proc.poll() is not None:
                        sleep(.2)
                        break
                    if line:
                        lane = line.decode().rstrip()
                        ml.stdout.append(lane)
                        LOG.debug(lane)
                        if lane.startswith('Elapsed time'):
                            break  # not parsable for iqtree2
                    else:
                        sleep(.2)
                        break

            if iface.pill2kill.is_set():
                sleep(.2)
                GObject.idle_add(self.stop_ML, errors, start)
                return True

            # check for errors
            for line in ml.stdout:
                if line.startswith('ERROR'):
                    errors.append(line)  # not parsable for iqtree2
            if errors:
                GObject.idle_add(self.stop_ML, errors, start)
                return True

            if mode == 'export':
                iface.text = 'building zip'
                LOG.debug(iface.text)
                sh = f"{'raxml' if ml.tool == 'raxml-ng' else 'iqtree'}_run.sh"
                with open(sh, 'w', newline='') as sf:
                    bash = ('''\n#!/bin/bash\n\n# Execute this script via 'bash '''
                            + ('raxml' if ml.tool == 'raxml-ng' else 'iqtree')
                            + '''_run.sh'\n
# get the number of CPUs available
cpus=$(nproc)\n
BLUE='\033[0;34m'
NC='\033[0m' # No Color
print_usage() {
  printf "Limit the number of threads/logical cores via ${BLUE}-f <number>${NC}. "\n}\n
cpu_limit=400\n
while getopts 'f:' flag; do
  case "${flag}" in
    f) cpu_limit="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;\n  esac\ndone\n
# find the minimum of the CPUs allowed and available
used=$(($cpu_limit<$cpus ? $cpu_limit : $cpus))\n
printf "${BLUE}$cpus${NC} CPUs available, use at most ${BLUE}$used${NC}.\nThis will proceed in a bit, interrupt with Ctrl+C\n"
print_usage
# make binary executable
chmod +x "'''
                            + ml.tool
                            + '''"\nsleep 10s''').strip()
                    sf.write(bash)
                    if ml.tool == 'raxml-ng':
                        format_args = {'raxml-ng': './raxml-ng',
                                       'msa': 'msa.fasta', 'seed': ml.ml_seed,
                                       'evo_model': Path(ml.evo_model).name,
                                       'model_modifier': ml.evo_modify,
                                       'chk_prefix': 'chk',
                                       'ml_prefix': 'ml',
                                       'bs_prefix': 'bs',
                                       'sp_prefix': 'sp',
                                       'cpus': '$used',
                                       'bootstraps': ml.bootstraps,
                                       'inferred_model': 'ml.raxml.bestModel',
                                       'ml_tree': 'ml.raxml.bestTree',
                                       'bs_trees_file': 'bs.raxml.bootstraps',
                                       'start_trees': ','.join(
                                           [a for a in ['rand{%d}' % ml.rand if ml.rand > 0 else None,
                                                        'pars{%d}' % ml.pars if ml.pars > 0 else None] if a])}
                        comments = ['\n\n# Check MSA\n', '\n# Find best ML tree\n',
                                    '\n# Compute bootstrap iterations\n',
                                    '\n# Calculate branch support\n']
                    elif ml.tool == 'iqtree2':
                        format_args = {'iqtree2': './iqtree2',
                                       'msa': 'msa.fasta', 'seed': ml.ml_seed,
                                       'evo_model': Path(ml.evo_model).name,
                                       'model_modifier': ml.evo_modify,
                                       'mf_prefix': 'mf',
                                       'ml_prefix': 'ml',
                                       'bs_prefix': 'bs',
                                       'uf_prefix': 'uf',
                                       'sp_prefix': 'sp',
                                       'cpus': '$used',
                                       'bootstraps': ml.bootstraps,
                                       'boot_trees': 'uf.ufboot' if ml.ultrafast else 'bs.boottrees',
                                       'start_trees': ml.rand + ml.pars}
                        comments = ['\n\n# Find best ML tree\n',
                                    '\n# Compute bootstrap iterations and FBP tree\n',
                                    '\n# Calculate TBE branch support\n']
                        if ml.infer_model:
                            # enable remote model inference
                            comments.insert(0, '\n\n# Find best model\n')
                            comments[1] = comments[1][1:]

                            format_args['evo_model'] = '$evo_model'
                            format_args['model_modifier'] = ''
                            format_args['BASH_REMATCH'] = [None, '{BASH_REMATCH[1]}']

                            calls.insert(0, '''
"{iqtree2}" -s "{msa}" -m MF -mtree -mset raxml -seed {seed} -nt AUTO -ntmax {cpus} -redo --prefix "mf"\n
echo "AB12PHYLO: iqtree2 ModelFinder finished"\n
# extract the suggested best model
line_regex="Best-fit model .*?$"\n
while read p; do
    # find the right line start
    if [[ $p =~ $line_regex ]]\n    then
        # revert the line to match the end; lazy modifiers don't work
        ledom_ove=$(echo $p | rev)
        model_regex="(\S+).*"
        if [[ $ledom_ove =~ $model_regex ]]
        then
            evo_model=$(echo ${BASH_REMATCH[1]} | rev)
            echo "AB12PHYLO: extracted model '$evo_model'"
        fi\n    fi\ndone < "mf.iqtree"'''.strip())
                    for comment, line in zip(comments, calls):
                        sf.write(comment)
                        sf.write(line.format(**format_args) + '\n')
                    sf.write('\n# Copy tree files\n')
                    if ml.tool == 'raxml-ng':
                        sf.write('cp sp.raxml.supportTBE tree_TBE.nwk\n')
                        sf.write('cp sp.raxml.supportFBP tree_FBP.nwk\n')
                    elif ml.tool == 'iqtree2':
                        sf.write(f'cp {"uf" if ml.ultrafast else "bs"}.contree tree_FBP.nwk\n')
                        sf.write(f'cp sp.suptree tree_TBE.nwk\n')
                    sf.write('\n')

                sleep(.05)
                GObject.idle_add(self.export_zip, sh, msa)
                return True

        if ml.in_shell and mode == 'run':  # same condition as at **1
            shell.chmod(shell.stat().st_mode | stat.S_IEXEC)
            self.hold()
            self.win.hide()
            sleep(.1)
            Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference running in background.',
                   '-i', str(repo.PATHS.icon_path)])
            call([shell])
            sleep(.1)
            self.win.show_all()
            self.release()

        iface.text = 'copy tree files'
        iface.i = iface.k - 1
        if ml.tool == 'raxml-ng':
            shutil.copy(prefix / 'sp.raxml.supportFBP', self.wd / repo.PATHS.fbp)
            shutil.copy(prefix / 'sp.raxml.supportTBE', self.wd / repo.PATHS.tbe)
        elif ml.tool == 'iqtree2':
            if ml.infer_model:
                # extract/parse model from mf.iqtree
                with open(f'{prefix / "mf"}.iqtree') as fh:
                    s = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
                    pos = s.find(b'Best-fit model')
                    if pos != -1:
                        model = s[pos + 10:pos + 190].decode().split('\n')[0].split(':')[-1].strip()
                        msg = f'Extracted model: {model}'
                        LOG.info(msg)
                        ml.evo_model = model.split('+')[0]
                        ml.evo_modify = model[len(ml.evo_model):]
                        if ml.evo_model not in data.evo_models.get_column(0):
                            data.evo_models.insert(len(data.evo_models) - 1, [ml.evo_model, None])
                        iface.evo_model.set_active_id(ml.evo_model)
                        iface.evo_modify.set_text(ml.evo_modify)
                    else:
                        errors.append(f'Could not extract model from\n{prefix / "mf"}.iqtree')
            else:
                shutil.copy(f'{prefix / ("uf" if ml.ultrafast else "bs")}.contree',
                            self.wd / repo.PATHS.fbp)
                shutil.copy(prefix / 'sp.suptree', self.wd / repo.PATHS.tbe)
        iface.text = 'idle'
        iface.frac = 1
        sleep(.1)
        GObject.idle_add(self.stop_ML, errors, start, msg)
        return True

    def stop_prep(self, errors):
        """Finish the ML prep thread"""
        iface = self.iface
        iface.thread.join()
        self.win.show_all()
        if errors:
            self.show_notification('Errors during ML prep', errors)
        iface.pill2kill.clear()
        self.refresh()
        return

    def stop_ML(self, errors, start, msg=False):
        """Finish the ML inference thread"""
        iface = self.iface
        iface.thread.join()
        self.update_ML(PAGE, self.data.ml)
        self.win.show_all()

        # change the button to its stop state
        im, la, tx = iface.tup
        im.set_from_icon_name('media-playback-start-symbolic', 4)
        la.set_text(tx)
        if errors:
            self.show_notification('Errors during ML inference', errors)
        elif iface.pill2kill.is_set():
            tx = 'stopped ML inference'
            self.show_notification(msg=tx, secs=2)
        elif time() - start > 120:
            Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference finished',
                   '-i', str(repo.PATHS.icon_path)])
            # notify = threading.Thread(target=_zenity, args=())
            # notify.start()
        else:
            tx = 'ML inference finished' if not msg else msg
            self.show_notification(tx)
        LOG.info(tx)
        iface.pill2kill.clear()
        self.refresh()
        if iface.run_after and not tx.startswith('stopped'):
            [do_func() for do_func in iface.run_after]
        self.set_changed(PAGE, False)
        return

    def _unset_modelfinder(self):
        self.iface.infer_model.set_active(False)
        self._load_tool_help(self.iface.ml_tool)

    def export_zip(self, sh, msa):
        """Finish the zip export thread"""
        data = self.data
        iface = self.iface
        ml = data.ml
        iface.thread.join()
        self.win.show_all()

        path = self.wd / f'{repo.toname(ml.tool)}_export.zip'
        dialog = Gtk.FileChooserDialog(title='export zip',
                                       parent=None, select_multiple=False,
                                       action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_folder(str(path.parent))
        dialog.set_current_name(path.name)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            try:
                p = Path(dialog.get_filename()).resolve()
                with ZipFile(p, 'w', ZIP_DEFLATED) as zf:
                    if ml.evo_model.endswith('user_model'):
                        zf.write(ml.evo_model, 'user_model')
                    # always pack the linux version of the used tool
                    zf.write(ml.binary if sys.platform == 'linux'
                             else ab12phylo_app_base.CFG.get(ml.tool + '-linux'), ml.tool)
                    zf.write(msa, 'msa.fasta')
                    zf.write(sh)
                Path(sh).unlink()
            except Exception as ex:
                LOG.error(ex)
        dialog.destroy()

        # change the button to its stop state
        im, la, tx = iface.tup
        im.set_from_icon_name('media-playback-start-symbolic', 4)
        la.set_text(tx)
