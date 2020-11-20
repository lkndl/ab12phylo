# 2020 Leo Kaindl

import logging
import threading
from pathlib import Path

import matplotlib.pyplot as plt

import gi
import numpy as np, sys
from time import sleep
from Bio import SeqIO
from argparse import Namespace
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib, GObject, GdkPixbuf

from GUI.gtk3 import commons, regex
from ab12phylo import filter

LOG = logging.getLogger(__name__)
plt.set_loglevel('warning')
PAGE = 2


# todo do not re-read if some were removed


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

    iface.view_qal.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)
    iface.view_qal.connect('size_allocate', set_dims, iface)
    iface.view_qal.connect('key_press_event', delete_and_ignore_rows, gui, PAGE, iface.view_qal.get_selection())
    iface.quality_refresh.connect('clicked', lambda *args: start_trim(gui))
    iface.quality_next.connect('clicked', commons.proceed, gui)
    iface.quality_back.connect('clicked', commons.step_back, gui)
    commons.bind_accelerator(gui.accelerators, iface.quality_next, '<Alt>Right')
    commons.bind_accelerator(gui.accelerators, iface.quality_back, '<Alt>Left')


def refresh(gui):
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
    start_trim(gui)


def start_trim(gui):
    """
    For non-empty gene set, start a background re-drawing thread
    and return the GUI to the main loop.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface
    # called after files are read
    if not data.genes or iface.running or iface.notebook.get_current_page() != PAGE:
        LOG.debug('abort re-draw')
        return

    if not data.seqdata:
        regex.start_read(gui, run_after=start_trim)
        return

    LOG.debug('start-up redraw')
    data.qal_model.clear()
    [child.destroy() for child in iface.qal_win.get_children()]
    sleep(.1)
    iface.thread = threading.Thread(target=do_trim, args=[gui])
    iface.running = True
    GObject.timeout_add(20, commons.update, iface, iface.plot_prog, PAGE)
    iface.thread.start()
    # GUI thread returns to main loop


def do_trim(gui):
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
    ignore_ids = iface.ignore_set if 'ignore_set' in iface else set()

    try:
        for record_id, gene in data.record_order:
            # skip records from other genes for the trimming preview
            if gene not in p.gene_roll or record_id in ignore_ids:
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
                row = commons.seqtoint(record)
            except AttributeError:
                # accept references anyway, but maybe skip no-phred ones
                if not data.metadata[gene][record_id]['is_ref'] and not p.accept_nophred:
                    continue
                has_qal, is_bad = False, False
                row = commons.seqtoint(record)
            except ValueError:
                has_qal, is_bad = True, True
                row = commons.seqtogray(record)
            rows.append(row)
            data.qal_model.append([record.id, has_qal, is_bad])

    except KeyError as ke:
        exit(ke)

    if not rows:
        LOG.warning('no sequence data remains')
    else:
        iface.txt = 'tabularize'
        max_len = max(map(len, rows))
        seq_array = np.array([row + commons.seqtoint(' ') * (max_len - len(row)) for row in rows])
        done += 1
        iface.frac = done / all_there_is_to_do

        iface.txt = 'plot'
        ratio = seq_array.shape[1] / seq_array.shape[0]
        masked = np.ma.masked_where(seq_array > commons.toint('N'), seq_array)
        f = Figure()  # figsize=(scale // 10 * ratio * 5, scale // 10), dpi=300)
        ax = f.add_subplot(111)
        mat = ax.matshow(masked, alpha=1, cmap=ListedColormap(commons.colors),
                         vmin=-.5, vmax=len(commons.colors) - .5, aspect='auto')
        ax.axis('off')
        f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        f.set_facecolor('none')

        # place on GTK
        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
        # canvas.mpl_connect('pick_event', onpick)
        # canvas.mpl_connect('button_press_event', onclick)

        done += 1
        iface.frac = done / all_there_is_to_do
        iface.txt = 'place + resize'
        data.width = seq_array.shape[1] * 4
        canvas.set_size_request(data.width, set_dims(iface.view_qal, None, iface)[1])
        # canvas.set_vexpand(False)
        try:
            iface.qal_win.add(canvas)
        except Gtk.Error as ex:
            LOG.error(ex)

        if gui.wd:
            if gui.wd:
                LOG.debug('saving plot')
                f.savefig(gui.wd / 'trim_preview.png', transparent=True,
                          dpi=600, bbox_inches='tight', pad_inches=0)

            with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
                LOG.debug('saving colorbar')
                fig = plt.figure(figsize=(4, 2))
                cax = fig.add_subplot(111)
                cbar = plt.colorbar(mat, ax=cax, ticks=range(len(commons.colors)), orientation='horizontal')
                cbar.ax.set_xticklabels(commons.bases)
                cax.remove()
                fig.savefig(gui.wd / 'colorbar.png', transparent=True,
                            bbox_inches='tight', pad_inches=0, dpi=600)
                del fig
            pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                str(gui.wd / 'colorbar.png'), 250, 100, preserve_aspect_ratio=True)
            iface.palplot.set_from_pixbuf(pb)

        iface.frac = 1

    sleep(.1)
    GObject.idle_add(stop_trim, gui)
    return True


def stop_trim(gui):
    iface = gui.interface
    iface.running = False
    iface.thread.join()
    gui.show_all()
    # link and resize scrollbar
    iface.qal_scroll.do_move_slider(iface.qal_scroll, Gtk.ScrollType.STEP_RIGHT)
    iface.plot_prog.set_text('idle')
    LOG.info('trim thread idle')
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
    """
    Edit a treeview cell in-place and save the result
    :param widget:
    :param data:
    :param iface:
    :return:
    """
    LOG.debug('editing')
    # filter for numbers only
    value = ''.join([c for c in widget.get_text() if c.isdigit()])
    widget.set_text(value if value else '')  # entering 0 is annoying


def delete_and_ignore_rows(widget, event, gui, page, selection):
    """
    Keep track of the rows that will not be written to the fasta and delete them from the treeview
    :param widget: 
    :param event: 
    :param gui: 
    :param page: 
    :param selection: 
    :return: 
    """
    data, iface = gui.data, gui.interface
    if Gdk.keyval_name(event.keyval) == 'Delete':
        model, iterator = selection.get_selected_rows()
        if 'ignore_set' not in iface:
            iface.ignore_set = set()
        for row in reversed(sorted(iterator)):
            iface.ignore_set.add(model[row[:]][0])
            model.remove(model.get_iter(row))

        commons.set_changed(gui, page, True)
        refresh(gui)


def parse(widget, event, gui):
    """
    Parse the content of a widget and if something changed cause a re-draw
    :param widget: The element to parse and inspect for changes
    :param event: sometimes passed by the signal
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface
    LOG.debug('parsing for re-draw')
    pre = None
    try:
        pre = iface.q_params.__getattribute__(widget.get_name())
    except KeyError:
        exit('%s not in q_params' % widget.get_name())
    # getting new value depends on widget type
    if widget == iface.gene_roll:
        now = widget.get_active_text()
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
    if pre == now:
        LOG.debug('no change, skip re-draw')
        return
    iface.q_params.__setattr__(widget.get_name(), now)
    if iface.q_params.trim_out > iface.q_params.trim_of:
        commons.show_message_dialog('cannot draw %d from %d' %
                                    (iface.q_params.trim_out, iface.q_params.trim_of))
        widget.set_text('0')
    else:
        commons.set_changed(gui, PAGE)
        start_trim(gui)


def trim_all(gui):
    """
    Trim all SeqRecords in project_dataset.seqdata to sequence strings and write collated .fasta files.
    Deleting SeqRecords will make the project file tiny again. If necessary, re-read trace files.
    Also place the png preview in the trim preview window.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.interface

    if not data.seqdata:
        LOG.debug('re-reading files')
        regex.start_read(gui, run_after=trim_all)
        return

    p = iface.q_params
    LOG.debug('writing collated .fasta files')
    for gene, genedata in data.seqdata.items():
        Path.mkdir(gui.wd / gene, exist_ok=True)

        # do actual trimming
        for _id in list(genedata.keys()):
            record = genedata.pop(_id)
            try:
                record = filter.trim_ends(record, p.min_phred, (p.trim_out, p.trim_of))
                record = filter.mark_bad_stretches(record, p.min_phred, p.bad_stretch)
            except ValueError:
                continue
            except AttributeError:
                pass
            genedata[_id] = record

        # write to file
        with open(str(gui.wd / gene / (gene + '.fasta')), 'w') as fasta:
            SeqIO.write(genedata.values(), fasta, 'fasta')

    # delete now bloaty data
    data.seqdata.clear()

    set_dims(iface.view_qal, None, iface)
    # place the png preview
    rectangle = iface.qal_win.get_allocated_size()[0]
    [child.destroy() for child in iface.qal_win.get_children()]
    pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(gui.wd / 'trim_preview.png'),
        width=data.width, height=rectangle.height,
        preserve_aspect_ratio=False)
    iface.qal_win.add(Gtk.Image.new_from_pixbuf(pixbuf))
    # link and resize scrollbar
    iface.qal_scroll.do_move_slider(iface.qal_scroll, Gtk.ScrollType.STEP_RIGHT)
    gui.show_all()


def set_dims(view_qal, event, iface):
    w, h = view_qal.get_allocated_width(), view_qal.get_allocated_height()
    iface.qal_spacer.set_size_request(w, -1)
    iface.qal_win.set_max_content_height(h)
    try:
        iface.qal_win.get_children()[0].set_size_request(w, h)
    except IndexError:
        pass
    return w, h


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
