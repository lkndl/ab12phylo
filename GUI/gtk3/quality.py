import logging
import queue, string
import pandas as pd
from pathlib import Path
from time import sleep
import threading

import gi
import re
import numpy as np
import seaborn as sns
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib, GObject

from GUI.gtk3 import commons
from ab12phylo import filter

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 2
VARIABLE = []


def init(gui):
    data, iface = gui.data, gui.interface
    iface.read_prog.set_visible(False)
    iface.gene_roll.set_model(Gtk.ListStore(str))
    iface.gene_roll.connect('changed', replace, data, iface)

    iface.accept_rev.set_active(iface.search_rev)
    iface.accept_nophred.set_active(True)

    iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=60, lower=0,
                                                  step_increment=1, page_increment=1))
    iface.min_phred.set_numeric(True)
    iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)

    iface.q_params = dict()

    for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
        wi = iface.__getattribute__(w_name)
        iface.q_params[w_name] = int(wi.get_text())
        wi.connect('changed', edit, data, iface)
        wi.connect('focus_out_event', parse, data, iface)
        wi.connect('activate', parse, None, data, iface)

    for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
        wi.connect('key-press-event', keypress, data, iface)

    for widget in [iface.accept_rev, iface.accept_nophred]:
        widget.connect('toggled', redraw, data, iface)

    iface.quality_next.connect('clicked', commons.proceed, gui)
    iface.quality_back.connect('clicked', commons.step_back, gui)

    reset(gui)


def reset(gui):
    if len(gui.data.rx_model) > 0:
        data, iface = gui.data, gui.interface

        iface.frac = 0
        iface.txt = ''

        iface.reader = threading.Thread(target=read, args=[(data, iface)])
        iface.running = True
        GObject.timeout_add(100, update, data, iface)
        iface.reader.start()
        print('when does this fire? EARLY')
        # GUI thread returns to main loop?


def update(data, iface):
    if iface.running:
        iface.notebook.get_children()[PAGE].set_sensitive(False)
        iface.read_prog.set_visible(True)
        iface.read_prog.set_fraction(iface.frac)
        iface.read_prog.set_text(iface.txt)
        return True
    else:
        iface.notebook.get_children()[PAGE].set_sensitive(True)
        iface.read_prog.set_visible(False)
        return False


def stop(iface, errors):
    iface.running = False
    iface.reader.join()
    iface.read_prog.set_text('idle')
    LOG.info('idle')
    if errors:
        commons.show_message_dialog('There were errors reading some files', errors)
    return


def read(args):
    data, iface = args
    data.csvs.clear()
    data.seqdata.clear()
    errors = list()
    do = len(data.rx_model) + len(data.wp_model)
    done = 0
    # read in wellsplates
    LOG.debug('reading wellsplates')
    iface.txt = 'reading plates ...'
    for row in data.wp_model:
        df = pd.read_csv(row[-1], header=None, engine='python')
        df.index = list(range(1, df.shape[0] + 1))
        df.columns = list(string.ascii_uppercase[0:df.shape[1]])
        box = row[0]
        if box in data.csvs:
            LOG.error('wellsplate %s already read in. overwrite with %s' % (box, row[-2]))
            commons.show_message_dialog(message='wellsplate %s already read in. overwrite with %s' % (box, row[-2]))
        data.csvs[box] = df
        done += 1
        iface.frac = done / do

    # read in trace files
    LOG.debug('reading traces')
    for row in data.rx_model:
        iface.txt = 'reading %s' % row[-2]
        file_path = row[-1]

        if file_path.endswith('.ab1'):
            # also check for phred scores. ABI traces also only contain a single record!
            sleep(0.16)
            try:
                record = SeqIO.read(file_path, 'abi')
            except UnicodeDecodeError:
                errors.append(row[-2])
            pass
        elif file_path.endswith('.fasta') or file_path.endswith('.fa') or file_path.endswith('.seq'):
            try:
                for record in SeqIO.parse(file_path, 'fasta'):
                    pass
            except UnicodeDecodeError:
                errors.append(row[-2])
        done += 1
        iface.frac = done / do
    # TODO start or call the initial redraw here
    GObject.idle_add(stop, iface, errors)
    return


def replace(widget, data, iface):
    LOG.debug('tabling ...')
    # transfer to GtkTreeView
    data.q_model.clear()

    redraw(None, data, iface)
    return


def redraw(widget, data, iface):
    LOG.debug('drawing ...')
    # get parameters from interface
    # also start this from re-sorting the treeview
    print(iface.q_params)


def delete_event(widget, event):
    return False


def keypress(widget, event, data, iface):
    key = Gdk.keyval_name(event.keyval)
    if key == 'Up':
        widget.set_text(str(1 + int(widget.get_text())))
        return True
    elif key == 'Down':
        widget.set_text(str(max(0, -1 + int(widget.get_text()))))
        return True


def edit(widget, data, iface):
    LOG.debug('editing ...')
    # filter for numbers only
    value = ''.join([c for c in widget.get_text() if c.isdigit()])
    widget.set_text(value if value else '0')


def parse(widget, event, data, iface):
    LOG.debug('parsing ...')
    # check if the numbers have changed
    pre = iface.q_params[widget.get_name()]
    now = int(widget.get_text())
    delete_event(widget, event)

    # cause redrawing if something changed
    if pre != now:
        iface.q_params[widget.get_name()] = now
        if iface.q_params['trim_out'] > iface.q_params['trim_of']:
            commons.show_message_dialog('cannot draw %d from %d' %
                                        (iface.q_params['trim_out'], iface.q_params['trim_of']))
            widget.set_text('0')
        else:
            redraw(None, data, iface)
