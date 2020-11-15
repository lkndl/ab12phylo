import logging
import queue, string
import pandas as pd
from pathlib import Path
from argparse import Namespace
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
KXLIN = [(.46, 1, .44, 1), (.16, .44, .8, 1), (1, .47, .66, 1), (.92, 1, .4, 1),
         (.84, .84, .84, .6), (1, 1, 1, 0), (1, 1, 1, 0), (1, 1, 1, 0), (1, 0, 0, 1)]


def init(gui):
    data, iface = gui.data, gui.interface
    iface.read_prog.set_visible(False)
    iface.gene_roll.set_model(data.gene_model)
    iface.gene_roll.connect('changed', parse, gui)

    iface.accept_rev.set_active(iface.search_rev)
    iface.accept_nophred.set_active(True)

    iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=60, lower=0,
                                                  step_increment=1, page_increment=1))
    iface.min_phred.set_numeric(True)
    iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)

    # this is kept up-to-date by the signals below
    iface.q_params = dict()

    for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
        wi = iface.__getattribute__(w_name)
        iface.q_params[w_name] = int(wi.get_text())
        wi.connect('changed', edit, data, iface)
        wi.connect('focus_out_event', parse, gui)
        wi.connect('activate', parse, None, gui)

    for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
        wi.connect('key-press-event', keypress, data, iface)

    for widget in [iface.accept_rev, iface.accept_nophred]:
        widget.connect('toggled', parse, gui)

    iface.qal_plot = Namespace()
    iface.qal_plot.sel = iface.view_dataset.get_selection()
    # TODO connect click inside plot with
    # iface.qal_plot.sel.select_iter(iter)

    # TODO init treeview
    # TODO sync treeview scrolling to plot scrolling

    iface.qal_plot.cmap = sns.color_palette(KXLIN, as_cmap=True)

    iface.quality_next.connect('clicked', commons.proceed, gui)
    iface.quality_back.connect('clicked', commons.step_back, gui)


def redraw(widget, gui):
    data, iface = gui.data, gui.interface
    # called after files are read
    # TODO start sth easier than this from re-sorting the treeview
    if not data.genes:
        return

    iface.bg_thread = threading.Thread(target=qplot, args=[gui])
    iface.running = True
    GObject.timeout_add(100, commons.update, iface, iface.plot_prog, PAGE)
    iface.bg_thread.start()
    # GUI thread returns to main loop


def qplot(gui):
    """
    Create a matrix representation of valid characters and plot as sns.heatmap
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface
    iface.frac = 0
    iface.txt = ''

    # parameters are up-to-date
    LOG.debug('re-draw with %s' % str(iface.q_params))

    # # create matrix from trace files
    # LOG.debug('create matrices')
    # iface.txt = 'sequence matrix'
    # data.seq_array = np.array([seq2ints(data.seqdata[_gene][_id].seq) for _id, _gene in data.record_ids])
    # LOG.debug('matrix shape: %s' % str(data.seq_array.shape))
    # done += 1
    # iface.frac = done / all_there_is_to_do
    # iface.txt = 'quality matrix'
    # data.qal_array = np.array([seq2qals(data.seqdata[_gene][_id], data.metadata[_gene][_id])
    #                            for _id, _gene in data.record_ids])
    # done += 1
    # iface.frac = done / all_there_is_to_do

    GObject.idle_add(stop, iface)
    return


def stop(iface):
    iface.running = False
    iface.bg_thread.join()
    iface.plot_prog.set_text('idle')
    LOG.info('idle')
    return


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


def parse(widget, event, gui):
    data, iface = gui.data, gui.interface
    LOG.debug('parsing ...')
    pre = iface.q_params[widget.get_name()]
    # getting new value depends on widget type
    if widget == iface.gene_roll:
        now = data.gene_model[event.get_active_iter()]
        now = data.genes if now == 'all' else set(now)
    elif widget in [iface.accept_rev, iface.accept_nophred]:
        now = widget.get_active()
    else:
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
            redraw(None, gui)
