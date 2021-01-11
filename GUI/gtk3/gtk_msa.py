# 2020 Leo Kaindl

import logging
import shutil
import subprocess
import threading
from argparse import Namespace
from pathlib import Path

import gi
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import GObject

from GUI.gtk3 import shared
from static import PATHS, TOOLS, toalgo

LOG = logging.getLogger(__name__)
PAGE = 3


def init(gui):
    data, iface = gui.data, gui.iface

    iface.msa_algo.set_entry_text_column(0)
    iface.msa_algo.set_id_column(0)
    iface.remote_algo.set_id_column(0)
    data.msa.cmd = dict()
    data.msa.remote_cmd = dict()
    data.msa.last_cmd = ''

    iface.msa_cmd.connect('focus_out_event', lambda widget, *args: data.msa.cmd.update(
        {toalgo(iface.msa_algo.get_active_text()): widget.get_buffer().props.text.strip()}))
    iface.remote_cmd.connect('focus_out_event', lambda widget, *args: data.msa.remote_cmd.update(
        {toalgo(iface.remote_algo.get_active_text()): widget.get_buffer().props.text.strip()}))

    iface.msa_algo.connect('changed', get_help, gui)
    iface.remote_algo.connect('changed', get_help, gui, True)
    iface.msa_import.connect('file-set', load_msa, gui)
    iface.msa_exe.connect('file-set', get_help, gui, True, True)

    # connect buttons
    iface.msa_build.connect('clicked', start_align, gui)
    shared.bind_accelerator(gui.accelerators, iface.msa_build, '<Enter>')
    iface.remote_build.connect('clicked', start_align, gui, True)
    shared.bind_accelerator(gui.accelerators, iface.remote_build, '<Enter>')


def get_help(widget, gui, remote=False, try_path=False):
    data, iface = gui.data, gui.iface

    if not data.genes:
        return
    if remote:
        data.msa.algo = toalgo(iface.remote_algo.get_active_text())
        client = TOOLS / 'MSA_clients' / (data.msa.algo + '.py')
        set_helpers(gui, 'python3 %s ' % client, iface.remote_help,
                    data.msa.remote_cmd, data.msa.algo, True, iface.remote_cmd)
    else:
        data.msa.algo = toalgo(iface.msa_algo.get_active_text())
        exe = shutil.which(data.msa.algo)
        exe = widget.get_active_text() if try_path else exe
        if exe:
            # get the --help output and save it in the lookup field on the right
            iface.msa_exe.set_filename(exe)
            set_helpers(gui, '%s --help; exit 0' % exe, iface.msa_help,
                        data.msa.cmd, data.msa.algo, False, iface.msa_cmd)
        else:
            # no executable found; unselect in path box
            iface.msa_exe.unselect_all()
            txt = data.msa.algo + ' was not found on your system $PATH. You can try ' \
                                  'manually specifying the path to the executable.'
            iface.msa_help.get_buffer().props.text = txt  # show this snarky line above
            iface.msa_cmd.get_buffer().props.text = ''  # no cmd suggestion


def set_helpers(gui, cmdline, help_view, help_dict, algo, remote, cmd_view):
    proc = subprocess.Popen(cmdline, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    txt = out.decode().strip() + '\n\n' + err.decode().strip()
    help_view.get_buffer().props.text = txt.strip()

    # get the suggested command and allow user modification in the left field
    txt = help_dict.get(algo, '')  # fetch saved
    if txt == '' or 'aligner' not in gui.iface:  # deleting all content will also get you back the original
        gui.iface.aligner, txt = shared.get_msa_build_cmd(algo, gui.wd, gui.data.genes, remote)  # fetch new
        help_dict[algo] = txt  # save

    cmd_view.get_buffer().props.text = help_dict[algo]  # show


def refresh_paths(gui):
    LOG.debug('refresh_paths')
    data, iface = gui.data, gui.iface
    old = str(iface.aligner.dir)
    iface.aligner.reset_paths(gui.wd, gui.wd / PATHS.msa)
    for algo, txt in data.msa.cmd.items():
        data.msa.cmd[algo] = txt.replace(old, str(gui.wd))
    # also replace it in the GtkTextView
    bf = iface.msa_cmd.get_buffer()
    bf.props.text = bf.props.text.replace(old, str(gui.wd))


def start_align(widget, gui, remote=False, run_after=None):
    """
    Starts an MSA building thread unless one of the following conditions is met:
    a) another thread is running -> abort + forbid proceeding.
    b) this function was called from the _Next button and an MSA already exists
    -> accept + allow proceeding.
    c) there were no changes registered for this page and no proceeding (not
    called by _Next) -> accept + show notification.

    :param widget: required for callback, ignored
    :param gui:
    :param remote: if the MSA should be constructed using the EMBL-EBI API
    :param run_after: the function to run afterwards; usually flip to next page
    :return:
    """
    data, iface = gui.data, gui.iface
    if iface.running:  # a)
        shared.show_notification(gui, 'Thread running')
        return
    elif data.msa.last_cmd == [data.msa.cmd, data.msa.remote_cmd][remote][data.msa.algo]:
        shared.show_notification(gui, 'MSA already generated, please proceed')
        return
    elif all([(gui.wd / gene / ('%s_raw_msa.fasta' % gene)).exists() for gene in data.genes]) \
            and run_after and not shared.get_errors(gui, PAGE):  # b)
        shared.set_changed(gui, PAGE, False)
        [do_func(gui) for do_func in run_after]
        return
    if 'aligner' not in iface:
        get_help(None, gui, remote)
    save_ui_state(gui)
    data.msa_lens.clear()
    iface.align_stack.props.sensitive = False
    iface.thread = threading.Thread(target=do_align, args=[gui, remote])
    iface.run_after = run_after
    iface.running = True
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    return
    # return to main loop


def do_align(gui, remote=False):
    data, iface = gui.data, gui.iface
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
                funcs[remote](gene, new_arg=arg % tuple([gene] * (4 - remote)))  # interpreting bool as int here
            except (FileNotFoundError, subprocess.CalledProcessError):
                refresh_paths(gui)
                # try again once more
                funcs[remote](gene, new_arg=arg % tuple([gene] * 4))
            # fetch MSA length
            for r in SeqIO.parse(gui.wd / gene / ('%s_raw_msa.fasta' % gene), 'fasta'):
                data.msa_lens.append(len(r))
                break
            iface.i += 1
        iface.frac = 1
        iface.text = 'idle'
    except (OSError, subprocess.CalledProcessError) as e:
        errors.append('%s at task %d (%s). invalid command?' % (type(e), iface.i, iface.text))
    except FileNotFoundError:
        errors.append('MSA/sequences file not found. Did you just save somewhere new?')
    except TypeError:
        errors.append('TODO replace string formatting with something smarter')
    GObject.idle_add(stop_align, gui, errors)


def stop_align(gui, errors):
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    iface.align_stack.props.sensitive = True
    gui.win.show_all()
    LOG.info('msa thread idle')
    shared.set_errors(gui, PAGE, bool(errors))
    shared.set_changed(gui, PAGE, False)
    if errors:
        shared.show_notification(gui, 'Errors during MSA building', errors)
        return
    if iface.run_after:
        [do_func(gui) for do_func in iface.run_after]
    else:
        shared.show_notification(gui, 'MSA building finished', items=None)
    return


def load_msa(widget, gui):
    data, iface = gui.data, gui.iface
    try:
        Path.mkdir(gui.wd / PATHS.import_msa.parent, exist_ok=True)
        shutil.copy(widget.get_filename(), gui.wd / PATHS.import_msa)
    except shutil.SameFileError:
        pass
    except Exception as ex:
        shared.show_notification(gui, str(ex))
        LOG.error(ex)
    shared.get_hashes(gui, PATHS.import_msa, PAGE)
    data.genes = ['import']
    data.gene_ids = {'import': {r.id for r in SeqIO.parse(gui.wd / PATHS.import_msa, 'fasta')}}
    iface.aligner, cmd = shared.get_msa_build_cmd(
        toalgo(gui.iface.msa_algo.get_active_text()), gui.wd, data.genes)
    LOG.debug('using imported MSA')
    save_ui_state(gui)


def refresh(gui):
    if gui.iface.align_stack.get_visible_child_name() == 'local':
        get_help(None, gui)
    elif gui.iface.align_stack.get_visible_child_name() == 'remote':
        get_help(None, gui, remote=True)


def save_ui_state(gui):
    ns, iface = gui.data.msa, gui.iface
    ns.stack_child_name = iface.align_stack.get_visible_child_name()
    ns.msa_algo_id = iface.msa_algo.get_active_id()
    ns.msa_exe_filename = iface.msa_exe.get_filename()
    ns.remote_algo_id = iface.remote_algo.get_active_id()
    ns.msa_import_filename = iface.msa_import.get_filename()


def reload_ui_state(gui):
    ns, iface = gui.data.msa, gui.iface
    iface.align_stack.set_visible_child_name(ns.stack_child_name)
    iface.msa_algo.set_active_id(ns.msa_algo_id)
    if ns.msa_exe_filename:
        iface.msa_exe.set_filename(ns.msa_exe_filename)
    iface.remote_algo.set_active_id(ns.remote_algo_id)
    if ns.msa_import_filename:
        iface.msa_import.set_filename(ns.msa_import_filename)

    iface.msa_cmd.get_buffer().props.text = ns.cmd.get(toalgo(ns.msa_algo_id), '')
    iface.remote_cmd.get_buffer().props.text = ns.remote_cmd.get(toalgo(ns.remote_algo_id), '')
