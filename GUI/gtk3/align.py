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
    data, iface = gui.data, gui.interface
    iface.plates = True

    iface.msa_algo.set_entry_text_column(0)
    iface.algos = Namespace()
    iface.algos.cmd = {algo: Gtk.TextBuffer() for algo in map(toalgo, algos)}

    iface.msa_cmd.connect('focus_out_event', lambda widget, *args: iface.algos.cmd.update(
        {toalgo(iface.msa_algo.get_active_text()): widget.get_buffer()}))

    iface.msa_handler = iface.msa_algo.connect('changed', get_help, gui)

    # connect buttons
    iface.align_next.connect('clicked', commons.proceed, gui)
    iface.align_back.connect('clicked', commons.step_back, gui)
    commons.bind_accelerator(gui.accelerators, iface.align_next, '<Alt>Right')
    commons.bind_accelerator(gui.accelerators, iface.align_back, '<Alt>Left')
    iface.msa_build.connect('clicked', start_align, gui)
    commons.bind_accelerator(gui.accelerators, iface.msa_build, '<Enter>')

    refresh(gui)


def get_help(widget, gui):
    data, iface = gui.data, gui.interface

    if not data.genes:
        return

    algo = toalgo(iface.msa_algo.get_active_text())

    exe = shutil.which(algo)
    if exe:
        # an executable was found
        iface.msa_exe.set_filename(exe)
        proc = subprocess.Popen(exe + ' --help; exit 0', shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        txt = out.decode().strip() + '\n\n' + err.decode().strip()
        iface.msa_help.get_buffer().props.text = txt

        buffer = iface.algos.cmd[algo]
        start, end = buffer.get_bounds()
        txt = buffer.get_text(start, end, True).strip()
        if txt == '':
            txt = get_cmd(algo, gui)
            buffer.props.text = txt

        iface.msa_cmd.props.buffer = iface.algos.cmd[algo]
    else:
        iface.msa_exe.unselect_all()
        txt = algo + ' was not found on your system $PATH. You can try manually' \
                     ' specifying the path to the executable.'
        iface.msa_help.get_buffer().props.text = txt
        iface.msa_cmd.props.buffer = Gtk.TextBuffer()


def get_cmd(algo, gui):
    data, iface = gui.data, gui.interface
    args = Namespace(**{
        'dir': gui.wd,
        'genes': data.genes,
        'msa_algo': algo,
        'user': commons.USER,
        'msa': gui.wd / 'msa.fasta',
        'sep': commons.SEP,
        'missing_samples': gui.wd / 'missing_samples.tsv'
    })
    iface.aligner = msa.msa_build(args, None, no_run=True)
    cmd = iface.aligner.build_local(data.agene(), True)
    return cmd


def start_align(widget, gui):
    data, iface = gui.data, gui.interface
    if 'aligner' not in iface:
        get_help(None, gui)
    iface.align_stack.props.sensitive = False
    iface.thread = threading.Thread(target=do_align, args=[gui])
    iface.run_after = None
    iface.running = True
    GObject.timeout_add(1, commons.update, iface, iface.align_prog, PAGE)
    iface.thread.start()
    # GUI thread returns to main loop


def do_align(gui):
    data, iface = gui.data, gui.interface
    exceptions = list()
    iface.frac = .05
    k = len(data.genes)
    i = 0
    for gene in data.genes:
        iface.txt = '%s [%d/%d]' % (gene, i + 1, k)
        try:
            iface.aligner.build_local(gene)
        except (OSError, subprocess.CalledProcessError) as e:
            exceptions.append(str(e))
        i += 1
        iface.frac = i / k
    GObject.idle_add(stop_align, gui, exceptions)


def stop_align(gui, errors):
    iface = gui.interface
    iface.running = False
    iface.thread.join()
    iface.align_stack.props.sensitive = True
    gui.show_all()
    iface.align_prog.props.text = 'idle'
    LOG.info('built MSAs')
    if errors:
        commons.show_message_dialog('Errors during MSA building', errors)
    return False


def refresh(gui):
    data, iface = gui.data, gui.interface
    if iface.align_stack.get_visible_child_name() == 'local':
        get_help(None, gui)
    # TODO elif
