# 2020 Leo Kaindl

import json
import logging
import re
import string
import sys
import shutil
import threading
import webbrowser
import shlex
import subprocess
from subprocess import run, Popen, PIPE
from pathlib import Path
from time import sleep

import gi
import pandas as pd
import requests, random
from Bio import SeqIO

import static

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import shared, gtk_qal
from ab12phylo import raxml
from static import PATHS, TOOLS

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 6


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    iface.evo_model.set_model(data.evo_models)
    iface.evo_model.get_child().connect('focus_out_event', _change_evo_model, data.ml)
    iface.evo_model.connect_after('changed', _load_model_file, iface.evo_modify, data.ml)
    iface.evo_model.set_active(0)

    iface.raxml_run.connect('clicked', start_ML, gui, 'raxml')

    for wi in [iface.bootstraps, iface.rand_trees, iface.pars_trees, iface.raxml_seed]:
        wi.connect('key-press-event', shared.edit_numerical_entry_up_down)
        wi.connect('changed', shared.edit_numerical_entry)

    # iface.raxml_run.connect('clicked', _load_raxml_help, gui, 'raxml')

    iface.raxml_seen = False


def _change_evo_model(entry, event_focus, ns):
    combo = entry.get_parent().get_parent()
    if combo.get_active_iter():
        ns.evo_model = entry.get_text()
    else:
        tx = entry.get_text()
        combo.get_model().append([tx, None])
        ns.evo_model = tx  # save in project dataset
        LOG.debug('entered custom evo model %s' % tx)


def _load_model_file(combo, evo_modify, ml):
    if combo.get_active_iter():
        if combo.get_active_id() != 'from file':
            evo_modify.props.sensitive = True
            ml.model_file = None
        else:
            dialog = Gtk.FileChooserDialog(title='select file with partition table',
                                           parent=None, select_multiple=False,
                                           action=Gtk.FileChooserAction.OPEN)
            dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                               Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = dialog.run()
            if response == Gtk.ResponseType.OK:
                tx = Path(dialog.get_filename()).resolve()
                tx, tx_path = tx.name, str(tx)
                combo.get_model().append([tx, tx_path])
                combo.set_active_id(tx)
                evo_modify.props.sensitive = False
                ml.model_file = tx_path
                LOG.debug('selected partitioned model file')
            dialog.destroy()


def _load_raxml_help(gui):
    data, iface = gui.data, gui.iface
    binary = shutil.which('raxml-ng')
    if not binary:
        binary = str(TOOLS / 'raxml-ng_v1.0.1_linux_x86_64' / 'raxml-ng')
    iface.raxml_exe.set_filename(binary)

    bf = iface.ml_help.get_buffer()

    with Popen(shlex.split('{} --help'.format(binary)),
               stdout=PIPE, shell=True) as proc:
        while True:
            line = proc.stdout.readline()
            if line:
                bf.insert(bf.get_end_iter(), line.decode(), -1)
                # bf.insert_at_cursor(line.decode(), -1)
            else:
                LOG.debug('got RAxML --help')
                break
    _draw_line(bf)


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree."""
    data, iface = gui.data, gui.iface

    if not iface.raxml_seen:
        _load_raxml_help(gui)
        iface.raxml_seen = True

    algo = static.toalgo(iface.ml_stack.get_visible_child_name())
    pass


def start_ML(widget, gui, mode):
    """Set-up the RAxML / IQ-Tree thread."""
    data, ml, iface = gui.data, gui.data.ml, gui.iface

    if not data.genes or iface.running:
        LOG.debug('abort ML')
        return

    if mode == 'raxml':
        ml.raxml_binary = iface.raxml_exe.get_filename()
        ml.evo_model = data.evo_models[iface.evo_model.get_active()][0]
        for w_name in ['evo_modify', 'bootstraps', 'rand_trees', 'pars_trees', 'raxml_seed']:
            ml.__setattr__(w_name, iface.__getattribute__(w_name).get_text())
            if w_name in ['bootstraps', 'rand_trees', 'pars_trees']:
                ml.__setattr__(w_name, int(ml.__getattribute__(w_name)))

        ml.raxml_seed = random.randint(0, max(1000, ml.bootstraps)) if ml.raxml_seed == '' else int(ml.raxml_seed)
        Path.mkdir(gui.wd / 'RAxML', exist_ok=True)
    else:
        raise NotImplementedError

    iface.thread = threading.Thread(target=do_ML, args=[gui, mode])
    iface.running = True
    GObject.timeout_add(1000, shared.update, iface, PAGE)
    iface.thread.start()
    # return to main loop


def do_ML(gui, mode):
    """Run the Gblocks thread."""
    data, ml, iface = gui.data, gui.data.ml, gui.iface
    bf = iface.ml_help.get_buffer()
    errors = list()
    iface.i = 0
    iface.k = ml.bootstraps + 2

    iface.text = 'checking MSA'
    output = ''
    with Popen(ml.raxml_binary + ' --msa ' + str(gui.wd / PATHS.msa) +
               ' --check --model ' + ml.evo_model + ml.evo_modify +
               ' --prefix ' + str(gui.wd / 'RAxML' / 'chk_'),
               stdout=PIPE, shell=True) as proc:
        while True:
            line = proc.stdout.readline()
            if line:
                line = line.decode()
                bf.insert(bf.get_end_iter(), line, -1)
                # bf.insert_at_cursor(line, -1)
                output += line
            else:
                break

    sleep(.1)
    GObject.idle_add(stop_ML, gui, errors)
    return True

    _draw_line(bf)
    iface.i += 1
    iface.text = 'start bootstrapping'
    output = ''
    with Popen([ml.raxml_binary, '--all', '--msa', gui.wd / PATHS.msa,
                '--model', ml.evo_model + ml.evo_modify,
                '--prefix', gui.wd / 'RAxML' / 'bs_',
                '--bs-trees autoMRE{%d}' % ml.bootstraps, '--seed', ml.raxml_seed,
                '--threads auto{16} --workers auto{8}', '--bs-metric fbp,tbe',
                '--tree rand{%d},pars{%d}' % (ml.rand_trees, ml.pars_trees)],
               stdout=PIPE, shell=True) as proc:
        while True:
            line = proc.stdout.readline()
            if line:
                line = line.decode()
                bf.insert_at_cursor(line, -1)
                output += line
            else:
                break

    _draw_line(bf)
    iface.text = 'idle'
    iface.frac = 1
    sleep(.1)
    GObject.idle_add(stop_ML, gui, errors)
    return True


def stop_ML(gui, errors):
    """Finish the ML inference thread"""
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    # TODO sensitivize a container?
    gui.win.show_all()
    LOG.info('ML thread idle')
    shared.set_errors(gui, PAGE, bool(errors))
    shared.set_changed(gui, PAGE, False)
    if errors:
        shared.show_notification(gui, 'Errors during ML inference', errors)
        return
    # if iface.run_after:
    #     [do_func(gui) for do_func in iface.run_after]
    return


def _draw_line(bf):
    bf.insert_markup(bf.get_end_iter(),
                     '<span foreground="#2374AF">'
                     '________________________________________________'
                     '________________________________\n</span>', -1)
