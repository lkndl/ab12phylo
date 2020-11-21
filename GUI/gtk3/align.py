# 2020 Leo Kaindl


import logging
import shutil
import subprocess
import threading
from argparse import Namespace

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import commons
from ab12phylo import msa

LOG = logging.getLogger(__name__)
PAGE = 3

algos = {'MAFFT': 'mafft', 'Clustal Omega': 'clustalo', 'MUSCLE': 'muscle', 'T-Coffee': 'tcoffee'}
toalgo = lambda c: algos[c]


def init(gui):
    data, iface = gui.data, gui.iface
    iface.plates = True

    iface.msa_algo.set_entry_text_column(0)
    iface.align = Namespace()
    iface.align.cmd = dict()
    iface.align.remote_cmd = dict()

    iface.msa_cmd.connect('focus_out_event', lambda widget, *args: iface.align.cmd.update(
        {toalgo(iface.msa_algo.get_active_text()): widget.get_buffer().props.text.strip()}))
    iface.remote_cmd.connect('focus_out_event', lambda widget, *args: iface.align.remote_cmd.update(
        {toalgo(iface.remote_algo.get_active_text()): widget.get_buffer().props.text.strip()}))

    iface.msa_algo.connect('changed', get_help, gui)
    iface.remote_algo.connect('changed', get_help, gui, True)
    iface.msa_import.connect('file-set', load_msa, gui)

    # connect buttons
    iface.msa_build.connect('clicked', start_align, gui)
    commons.bind_accelerator(gui.accelerators, iface.msa_build, '<Enter>')
    iface.remote_build.connect('clicked', start_align, gui, True)
    commons.bind_accelerator(gui.accelerators, iface.remote_build, '<Enter>')


def get_help(widget, gui, remote=False):
    data, iface = gui.data, gui.iface

    if not data.genes:
        return

    algo = toalgo(iface.msa_algo.get_active_text())
    if remote:
        client = commons.TOOLS / 'MSA_clients' / (algo + '.py')
        set_helpers(gui, 'python3 %s ' % client, iface.remote_help,
                    iface.align.remote_cmd, algo, True, iface.remote_cmd)
    else:
        exe = shutil.which(algo)
        if exe:
            # get the --help output and save it in the lookup field on the right
            iface.msa_exe.set_filename(exe)
            set_helpers(gui, '%s --help; exit 0' % exe, iface.msa_help,
                        iface.align.cmd, algo, False, iface.msa_cmd)
        else:
            # no executable found; unselect in path box
            iface.msa_exe.unselect_all()
            txt = algo + ' was not found on your system $PATH. You can try manually' \
                         ' specifying the path to the executable.'
            iface.msa_help.get_buffer().props.text = txt  # show this snarky line
            iface.msa_cmd.get_buffer().props.text = ''  # no cmd suggestion


def set_helpers(gui, cmdline, help_view, help_dict, algo, remote, cmd_view):
    proc = subprocess.Popen(cmdline, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    txt = out.decode().strip() + '\n\n' + err.decode().strip()
    help_view.get_buffer().props.text = txt

    # get the suggested command and allow user modification in the left field
    txt = help_dict.get(algo, '')  # fetch saved
    if txt == '':  # deleting all content will also get you back the original
        txt = get_cmd(algo, gui, remote)  # fetch new
        help_dict[algo] = txt  # save

    cmd_view.get_buffer().props.text = help_dict[algo]  # show


def get_cmd(algo, gui, remote=False):
    data, iface = gui.data, gui.iface
    args = Namespace(**{
        'dir': gui.wd,
        'genes': list(data.genes),
        'msa_algo': algo,
        'user': commons.USER,
        'msa': gui.wd / commons.RAW_MSA,
        'sep': commons.SEP,
        'missing_samples': gui.wd / 'missing_samples.tsv'
    })
    iface.aligner = msa.msa_build(args, None, no_run=True)
    if remote:
        cmd = iface.aligner.build_remote(data.agene(), no_run=True)
    else:
        cmd = iface.aligner.build_local(data.agene(), no_run=True)
    return cmd


def start_align(widget, gui, remote=False):
    data, iface = gui.data, gui.iface
    if 'aligner' not in iface:
        get_help(None, gui, remote)
    iface.align_stack.props.sensitive = False
    iface.thread = threading.Thread(target=do_align, args=[gui, remote])
    iface.run_after = None
    iface.running = True
    GObject.timeout_add(1000, commons.update, iface, PAGE)
    iface.thread.start()
    # return to main loop


def do_align(gui, remote=False):
    data, iface = gui.data, gui.iface
    exceptions = list()
    iface.frac = .05
    k = len(data.genes) + 2
    i = 0
    for gene in data.genes:
        iface.txt = 'aligning %s [%d/%d]' % (gene, i + 1, k - 1)
        try:
            if remote:
                iface.aligner.build_remote(gene)
            else:
                iface.aligner.build_local(gene)
        except (OSError, subprocess.CalledProcessError) as e:
            exceptions.append(str(e))
        i += 1
        iface.frac = i / k
    LOG.info('built MSAs')
    iface.txt = 'concatenating MSAs'
    iface.aligner.concat_msa(raw=True)
    i += 1
    iface.frac = i / k
    iface.txt = 'comparing SHA256 hashes'
    compare_hashes(gui)
    GObject.idle_add(stop_align, gui, exceptions)


def stop_align(gui, errors):
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    iface.align_stack.props.sensitive = True
    gui.win.show_all()
    iface.prog_bar.props.text = 'idle'
    LOG.info('align thread idle')
    if errors:
        commons.show_notification(gui, 'Errors during MSA building', errors)
    return False


def load_msa(widget, gui):
    try:
        shutil.copy(widget.get_filename(), gui.wd / commons.RAW_MSA)
        compare_hashes(gui)
    except shutil.SameFileError:
        pass
    except Exception as ex:
        commons.show_notification(gui, str(ex))
        LOG.error(ex)


def refresh(gui):
    data, iface = gui.data, gui.iface
    if iface.align_stack.get_visible_child_name() == 'local':
        get_help(None, gui)
    elif iface.align_stack.get_visible_child_name() == 'remote':
        get_help(None, gui, remote=True)


def compare_hashes(gui):
    msa_hash = commons.file_hash(gui.wd / commons.RAW_MSA)
    if msa_hash != gui.data.msa_hash:
        commons.set_changed(gui, PAGE)
        gui.data.msa_hash = msa_hash
