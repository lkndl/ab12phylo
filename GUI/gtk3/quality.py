import logging
from pathlib import Path
import threading
import sys

import gi
import numpy as np
import seaborn as sns
from time import sleep
from argparse import Namespace
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.figure import Figure

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib, GObject

from GUI.gtk3 import commons
from ab12phylo import filter

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 2
KXLIN = [('A', (0.46, 1, 0.44, 1)),
         ('C', (0.16, 0.44, 0.8, 1)),
         ('G', (1, 0.47, 0.66, 1)),
         ('T', (0.92, 1, 0.4, 1)),
         ('N', (0.84, 0.84, 0.84, 0.6)),
         ('-', (1, 1, 1, 0)),
         ('[ ]', (1, 1, 1, 0)),
         ('sep', (1, 1, 1, 0)),
         ('unknown', (1, 0, 0, 1))]
CMAP = sns.color_palette([c[1] for c in KXLIN], as_cmap=True)


def init(gui):
    data, iface = gui.data, gui.interface
    iface.read_prog.set_visible(False)
    iface.gene_roll.set_entry_text_column(0)
    iface.gene_handler = iface.gene_roll.connect('changed', parse, None, gui)
    iface.rev_handler = iface.accept_rev.connect('toggled', parse, None, gui)
    iface.accept_nophred.set_active(True)
    iface.accept_nophred.connect('toggled', parse, None, gui)

    iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=60, lower=0,
                                                  step_increment=1, page_increment=1))
    iface.min_phred.set_numeric(True)
    iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)
    iface.qal_scroll.connect('change_value', lambda *args: iface.
                             qal_win.set_hadjustment(args[0].get_adjustment()))

    # this is kept up-to-date by the signals below
    iface.q_params = Namespace(gene_roll='all')

    for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
        wi = iface.__getattribute__(w_name)
        iface.q_params.__setattr__(w_name, int(wi.get_text()))
        wi.connect('changed', edit, data, iface)
        wi.connect('focus_out_event', parse, gui)
        wi.connect('activate', parse, None, gui)

    for w_name in ['accept_rev', 'accept_nophred']:
        iface.q_params.__setattr__(w_name, iface.__getattribute__(w_name).get_active())

    for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
        wi.connect('key-press-event', keypress, data, iface)

    # init row annotation
    iface.view_qal.set_model(data.qal_model)
    iface.view_qal.set_headers_visible(False)
    iface.view_qal.append_column(Gtk.TreeViewColumn(
        title='id', cell_renderer=Gtk.CellRendererText(), text=0))
    iface.view_qal.append_column(Gtk.TreeViewColumn(
        title='has phreds', active=1,
        cell_renderer=Gtk.CellRendererToggle(radio=False)))
    iface.view_qal.connect('size_allocate', set_dims, iface)

    iface.quality_next.connect('clicked', commons.proceed, gui)
    iface.quality_back.connect('clicked', commons.step_back, gui)


def reset(gui):
    data, iface = gui.data, gui.interface
    with iface.gene_roll.handler_block(iface.gene_handler):
        iface.gene_roll.remove_all()

        [iface.gene_roll.append_text(gene) for gene in data.genes]
        if len(data.genes) > 1:
            iface.gene_roll.insert_text(0, 'all')
        iface.gene_roll.set_active(0)
    with iface.accept_rev.handler_block(iface.rev_handler):
        iface.accept_rev.set_active(iface.reverse_rx_chk.get_active())
        iface.accept_rev.set_sensitive(iface.reverse_rx_chk.get_active())

    redraw(gui)


def redraw(gui):
    """
    For non-empty gene set, start a background re-drawing thread
    and return the GUI to the main loop.
    :param gui:
    :return:
    """
    LOG.debug('start-up redraw')
    data, iface = gui.data, gui.interface
    # called after files are read
    # TODO start sth easier than this from re-sorting the treeview
    # TODO palplot in upper half

    # (gtk_main.py:5417): Gtk-CRITICAL **: 22:28:20.481: gtk_container_forall: assertion 'GTK_IS_CONTAINER (container)' failed
    if not data.genes or iface.running or iface.notebook.get_current_page() != PAGE:
        LOG.debug('abort re-draw')
        return

    data.qal_model.clear()
    [child.destroy() for child in iface.qal_win.get_children()]
    sleep(.1)
    iface.thread = threading.Thread(target=qplot, args=[gui])
    iface.running = True
    GObject.timeout_add(20, commons.update, iface, iface.plot_prog, PAGE)
    iface.thread.start()
    # GUI thread returns to main loop


def qplot(gui):
    """
    Iterate over records and trim, create a matrix representation
    of valid characters and plot as seaborn heatmap.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface
    # parameters are up-to-date
    LOG.debug('re-draw with %s' % str(iface.q_params))
    iface.frac = 0
    iface.txt = 'creating matrix'
    rows = list()
    p = iface.q_params
    p.gene_roll = data.genes if p.gene_roll == 'all' else p.gene_roll
    done = 0
    all_there_is_to_do = sum([len(data.seqdata[gene]) for gene in p.gene_roll]) + 3

    try:
        for record_id, gene in data.record_order:
            # skip records from other genes for the trimming preview
            if gene not in p.gene_roll:
                continue

            done += 1
            iface.frac = done / all_there_is_to_do
            # maybe skip reversed seqs
            if not p.accept_rev and data.metadata[gene][record_id]['is_rev']:
                continue

            record = data.seqdata[gene][record_id]
            try:
                record = filter.trim_ends(record, p.min_phred, (p.trim_out, p.trim_of), trim_preview=True)
                record = filter.mark_bad_stretches(record, p.min_phred, p.bad_stretch)
                has_qal, is_bad = True, False
                row = filter.seq2ints(record)
            except AttributeError:
                # accept references anyway, but maybe skip no-phred ones
                if not data.metadata[gene][record_id]['is_ref'] and not p.accept_nophred:
                    continue
                has_qal, is_bad = False, False
                row = filter.seq2ints(record)
            except ValueError:
                has_qal, is_bad = True, True
                row = filter.seq2gray(record)
            rows.append(row)
            data.qal_model.append([record.id, has_qal, is_bad])

    except KeyError as ke:
        exit(ke)

    if not rows:
        exit('no sequence data remains')
    else:
        set_dims(iface.view_qal, None, iface)

        iface.txt = 'tabularize'
        max_len = max(map(len, rows))
        seq_array = np.array([row + [5] * (max_len - len(row)) for row in rows])  # MARK 5 is the gap character
        done += 1
        iface.frac = done / all_there_is_to_do

        iface.txt = 'plot'
        sns.set(style='white')
        sns.despine(offset=0, trim=True)
        ratio = seq_array.shape[1] / seq_array.shape[0]
        figure = Figure(dpi=72)

        sns.heatmap(seq_array, ax=figure.add_subplot(111), cmap=CMAP,  # DO NOT USE annot=True,
                    yticklabels=False, xticklabels=False,  # DO NOT USE square=True
                    vmin=-.5, vmax=len(KXLIN) - .5,  # adjust the color map to the character range
                    linewidth=0, cbar=False)
        figure.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        canvas = FigureCanvas(figure)  # a Gtk.DrawingArea
        # canvas.mpl_connect('pick_event', onpick)
        # canvas.mpl_connect('button_press_event', onclick)

        done += 1
        iface.frac = done / all_there_is_to_do
        iface.txt = 'place + resize'
        canvas.set_size_request(width=max(4 * seq_array.shape[1], p.height * ratio // 2),
                                height=p.height)
        try:
            iface.qal_win.add(canvas)
        except Gtk.Error as ex:
            LOG.error(ex)
        iface.frac = .99
        # plt.ion()

    sleep(.1)
    GObject.idle_add(stop, gui)
    return True


def stop(gui):
    iface = gui.interface
    iface.running = False
    iface.thread.join()
    # del iface.thread
    gui.show_all()
    set_dims(iface.view_qal, None, iface)
    # link and resize scrollbar
    iface.qal_scroll.do_move_slider(iface.qal_scroll, Gtk.ScrollType.STEP_RIGHT)
    iface.plot_prog.set_text('idle')
    LOG.info('idle')
    return False


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
    widget.set_text(value if value else '')  # entering 0 is annoying


def parse(widget, event, gui):
    """
    Parse the content of a widget and if something changed cause a re-draw
    :param widget: The element to parse and inspect for changes
    :param event: sometimes passed by the signal
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface
    print(event, file=sys.stderr)
    print(widget, file=sys.stderr)
    LOG.debug('parsing for re-draw')
    pre = None
    try:
        pre = iface.q_params.__getattribute__(widget.get_name())
    except KeyError:
        exit('%s not in q_params' % widget.get_name())
    # getting new value depends on widget type
    if widget == iface.gene_roll:
        now = widget.get_active_text()
        # if now == {}:
        #     return
        now = data.genes if now == 'all' else {now}
    elif widget in [iface.accept_rev, iface.accept_nophred]:
        now = widget.get_active()
    else:
        try:
            now = int(widget.get_text())
        except ValueError:
            now = 0
    delete_event(widget, event)

    # cause re-drawing if something changed
    if pre != now:
        iface.q_params.__setattr__(widget.get_name(), now)
        if iface.q_params.trim_out > iface.q_params.trim_of:
            commons.show_message_dialog('cannot draw %d from %d' %
                                        (iface.q_params.trim_out, iface.q_params.trim_of))
            widget.set_text('0')
        else:
            redraw(gui)
    else:
        LOG.debug('no change, abort re-draw')
    return True


def trim_all(gui):
    """
    Really trim all sequences and write to new data structure.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface


def set_dims(widget, rect, iface):
    iface.q_params.height = min(20, widget.get_allocated_height())
    iface.qal_spacer.set_size_request(widget.get_allocated_width(), -1)
    return True


def onpick(event):
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = tuple(zip(xdata[ind], ydata[ind]))
    print('onpick points:', points)
    print('This would be useful after the next big restructuring', file=sys.stderr)
    return True


def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
    print('This would be useful after the next big restructuring', file=sys.stderr)
    return True
