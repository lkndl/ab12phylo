# 2020 Leo Kaindl

import logging
import threading
from argparse import Namespace
from pathlib import Path
from time import sleep

import gi
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.colorbar import ColorbarBase

import static

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GObject

from static import PATHS
from GUI.gtk3 import shared, gtk_rgx

from ab12phylo.filter import trim_ends, mark_bad_stretches

LOG = logging.getLogger(__name__)
plt.set_loglevel('warning')
PAGE = 2

"""The page for ABI trace trimming. Very similar to the MSA trimming page gtk_gbl.py"""


def init(gui):
    data, iface = gui.data, gui.iface

    iface.accept_rev.set_active(data.accept_rev)
    iface.rev_handler = iface.accept_rev.connect('toggled', parse, None, gui)
    iface.accept_nophred.set_active(True)
    iface.accept_nophred.connect('toggled', parse, None, gui)

    iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=61, lower=0,
                                                  step_increment=1, page_increment=1))
    iface.min_phred.set_numeric(True)
    iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)

    # this is kept up-to-date by the signals below
    iface.qal = Namespace(gene_roll='all')

    for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
        wi = iface.__getattribute__(w_name)
        iface.qal.__setattr__(w_name, int(wi.get_text()))
        wi.connect('changed', _edit_numerical_entry)
        wi.connect('focus_out_event', parse, gui)
    iface.min_phred.disconnect_by_func(parse)  # MARK leave this line alone

    for w_name in ['accept_rev', 'accept_nophred']:
        iface.qal.__setattr__(w_name, iface.__getattribute__(w_name).get_active())

    for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
        wi.connect('key-press-event', _edit_numerical_entry_up_down)

    # init row annotation
    iface.view_qal.set_model(data.qal_model)
    iface.view_qal.set_headers_visible(False)
    iface.view_qal.append_column(Gtk.TreeViewColumn(
        title='id', cell_renderer=Gtk.CellRendererText(), text=0, underline=2, strikethrough=3))
    # crt = Gtk.CellRendererToggle(radio=False)
    # crt.props.indicator_size = 13
    # iface.view_qal.append_column(Gtk.TreeViewColumn(
    #     title='no phreds', active=1, cell_renderer=crt))

    sel = iface.view_qal.get_selection()
    sel.set_mode(Gtk.SelectionMode.MULTIPLE)
    sel.connect('changed', shared.keep_visible,
                iface.parallel_qal.props.vadjustment.props, iface.qal)
    iface.view_qal.connect('check-resize', shared.get_height_resize,
                           iface.qal_spacer, [iface.qal_win])
    # in-preview deletion
    iface.view_qal.connect('key_press_event', shared.delete_and_ignore_rows,
                           gui, PAGE, sel, iface.qal)
    iface.qal_eventbox.connect_after('button_press_event', shared.select_seqs, PAGE, iface.zoom,
                                     iface.view_qal, iface.qal)  # in-preview selection
    iface.qal_win.connect('scroll-event', shared.xy_scale, gui, PAGE)  # zooming


def refresh(gui):
    data, iface = gui.data, gui.iface
    if not (gui.wd / PATHS.preview).exists() or 0 in data.qal_shape:
        start_trim(gui)
        return

    with iface.accept_rev.handler_block(iface.rev_handler):
        iface.accept_rev.set_active(iface.reverse_rx_chk.get_active() and iface.accept_rev.get_active())
        iface.accept_rev.set_sensitive(iface.reverse_rx_chk.get_active())

    # place the png preview
    shared.load_image(iface.zoom, PAGE, iface.qal_eventbox, gui.wd / PATHS.preview,
                      data.qal_shape[0] * shared.get_hadj(iface), data.qal_shape[1])
    shared.load_colorbar(iface.palplot, gui.wd)
    gui.win.show_all()


def parse(widget, event, gui):
    """
    Parse the content of a widget. Hitting enter in an entry still calls this, which is good.
    :param widget: The element to parse and inspect for changes
    :param event: passed by some signals -> focus_out_event
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.iface
    LOG.debug('parsing for re-draw')
    pre = None
    try:
        pre = iface.qal.__getattribute__(widget.get_name())
    except KeyError:
        exit('%s not in q_params' % widget.get_name())
    # getting new value depends on widget type
    if widget == iface.gene_roll:
        now = widget.get_active_text()
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
        LOG.debug('no change at parsing')
        return
    iface.qal.__setattr__(widget.get_name(), now)
    if iface.qal.trim_out > iface.qal.trim_of:
        shared.show_notification(gui, 'cannot draw %d from %d' %
                                 (iface.qal.trim_out, iface.qal.trim_of))
        widget.set_text('0')
    else:
        shared.set_changed(gui, PAGE)


def start_trim(gui):
    """
    For non-empty gene set, start a background re-drawing thread
    and return the GUI to the main loop.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.iface
    # called after files are read
    if not data.genes or iface.running or iface.notebook.get_current_page() != PAGE:
        LOG.debug('abort re-draw')
        return

    if not data.seqdata:
        gtk_rgx.start_read(gui, run_after=[start_trim])
        return

    data.accept_rev = iface.accept_rev.get_active()
    parse(iface.min_phred, None, gui)  # annoying SpinButton

    LOG.debug('start-up redraw')
    data.qal_model.clear()
    sleep(.1)
    iface.thread = threading.Thread(target=do_trim, args=[gui])
    iface.running = True
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    # return to main loop


def do_trim(gui):
    """
    Iterate over records and trim, create a matrix representation
    of valid characters and plot as seaborn heatmap.
    :param gui:
    :return:
    """
    data, iface = gui.data, gui.iface
    # parameters are up-to-date
    LOG.debug('re-draw with %s' % str(iface.qal))
    p = iface.qal
    rows = list()
    data.gene_for_preview = p.gene_roll
    genes_for_preview = data.genes if p.gene_roll == 'all' else {p.gene_roll}
    ignore_ids = iface.qal.ignore_ids if 'ignore_ids' in iface.qal \
        else {g: dict() for g in data.genes}  # regulated at delete_and_ignore_rows
    iface.text = 'creating matrix'
    LOG.debug(iface.text)
    iface.i = 0
    iface.k = sum([len(data.seqdata[gene]) for gene in genes_for_preview]) + 4  # number of all records + extra
    shared_ids = set.intersection(*data.gene_ids.values())

    try:
        for record_id, gene in data.record_order:
            # skip records from other genes for the trimming preview
            if gene not in genes_for_preview or \
                    gene in ignore_ids and record_id in ignore_ids[gene]:
                continue

            iface.i += 1
            # maybe skip reversed seqs
            if not p.accept_rev and data.metadata[gene][record_id]['is_rev']:
                continue

            record = data.seqdata[gene][record_id]
            try:
                record = trim_ends(record, p.min_phred, (p.trim_out, p.trim_of), trim_preview=True)
                record = mark_bad_stretches(record, p.min_phred, p.bad_stretch)
                has_qal, is_bad = True, False
                row = static.seqtoint(record)
            except AttributeError:
                # accept references anyway, but maybe skip no-phred ones
                is_ref = 'accession' in data.metadata[gene][record_id]
                if not is_ref and not p.accept_nophred:
                    continue
                has_qal, is_bad = is_ref, False
                row = static.seqtoint(record)
            except ValueError:
                has_qal, is_bad = True, True
                row = static.seqtogray(record)
            rows.append(row)
            underline = not has_qal if record.id in shared_ids else 4
            data.qal_model.append([record.id, gene, underline, is_bad])

    except KeyError as ke:
        LOG.error(ke)

    if not rows:
        LOG.warning('no sequence data remains')
        sleep(.1)
        GObject.idle_add(stop_trim, gui)
        return True

    iface.text = 'tabularize'
    LOG.debug(iface.text)
    max_len = max(map(len, rows))
    array = np.array([row + static.seqtoint(' ') * (max_len - len(row)) for row in rows])
    # make gaps transparent
    array = np.ma.masked_where(array > static.toint('else'), array)
    iface.i += 1
    iface.text = 'plot'
    LOG.debug(iface.text)

    # adjust maximum size
    scale = 6
    while max(array.shape) * scale > 2 ** 14:
        scale -= 1
    LOG.debug('scaling qal with %d' % scale)

    f = Figure()
    f.set_facecolor('none')
    f.set_figheight(array.shape[0] / static.DPI)
    f.set_figwidth(array.shape[1] / static.DPI)
    ax = f.add_subplot(111)
    ax.axis('off')
    f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    mat = ax.matshow(array, alpha=1, cmap=ListedColormap(static.colors),
                     vmin=-.5, vmax=len(static.colors) - .5, aspect='auto')

    iface.i += 1
    iface.text = 'save PNG'
    LOG.debug(iface.text)
    Path.mkdir(gui.wd / PATHS.preview.parent, exist_ok=True)
    f.savefig(gui.wd / PATHS.preview, transparent=True,
              dpi=scale * static.DPI, bbox_inches='tight', pad_inches=0.00001)
    plt.close(f)

    data.qal_shape[0] = array.shape[1]
    data.qal_shape[1] = shared.get_height_resize(iface.view_qal, None, iface.qal_spacer, [iface.qal_win])

    if iface.rasterize.props.active:
        iface.text = 'place PNG'
        LOG.debug(iface.text)
        sleep(.05)
        shared.load_image(iface.zoom, PAGE, iface.qal_eventbox, gui.wd / PATHS.preview,
                          data.qal_shape[0] * shared.get_hadj(iface), data.qal_shape[1])
    else:
        iface.text = 'place vector'
        LOG.debug(iface.text)
        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
        canvas.set_size_request(data.qal_shape[0] * shared.get_hadj(iface), data.qal_shape[1])
        try:
            iface.qal_eventbox.get_child().destroy()
            iface.qal_eventbox.add(canvas)
            # iface.qal_tools.attach(NavigationToolbar(canvas, gui.win), 1, 1, 1, 1)
        except Gtk.Error as ex:
            LOG.error(ex)
    iface.i += 1

    with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
        iface.text = 'colorbar'
        LOG.debug(iface.text)
        fig = plt.figure(figsize=(4, .5))
        cax = fig.add_subplot(111)
        i = static.NUCLEOTIDES.index('-')
        cbar = ColorbarBase(ax=cax, cmap=ListedColormap(static.colors[:i]),
                            ticks=[(j + .5) / i for j in range(i)], orientation='horizontal')
        cbar.ax.set_xticklabels(static.NUCLEOTIDES[:i])
        fig.savefig(gui.wd / PATHS.cbar, transparent=True,
                    bbox_inches='tight', pad_inches=0, dpi=600)
        plt.close(fig)
        del fig, cbar
    shared.load_colorbar(iface.palplot, gui.wd)

    iface.text = 'idle'
    iface.frac = 1
    sleep(.1)
    GObject.idle_add(stop_trim, gui)
    return True


def stop_trim(gui):
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    gui.win.show_all()
    iface.prog_bar.set_text('idle')
    LOG.info('qal thread idle')
    return False


def delete_event(widget, event):
    return False


def _edit_numerical_entry_up_down(widget, event):
    key = Gdk.keyval_name(event.keyval)
    if key == 'Up':
        widget.set_text(str(1 + int(widget.get_text())))
        return True
    elif key == 'Down':
        widget.set_text(str(max(0, -1 + int(widget.get_text()))))
        return True


def _edit_numerical_entry(widget):
    LOG.debug('editing')
    # filter for numbers only
    value = ''.join([c for c in widget.get_text() if c.isdigit()])
    widget.set_text(value if value else '')  # entering 0 is annoying


def trim_all(gui, run_after=None):
    """
    Trim all SeqRecords in project_dataset.seqdata to sequence strings
    and write two collated .fasta files per gene, one only containing
    the records shared across all genes. If necessary, trace files are
    re-read before and will be deleted from memory afterwards for a
    smaller project file. Plot and place png previews of trimming result.
    :param gui:
    :param run_after: [functions(gui)] when finished, usually flip page.
    :return:
    """
    data, iface = gui.data, gui.iface

    if not data.seqdata:
        LOG.debug('re-reading files')
        gtk_rgx.start_read(gui, run_after=[trim_all])
        return

    p = iface.qal
    ignore_ids = iface.qal.ignore_ids if 'ignore_ids' in iface.qal \
        else {g: dict() for g in data.genes}  # regulated at delete_and_ignore_rows
    LOG.debug('trim and filter all sequences')
    for gene, genedata in data.seqdata.items():
        Path.mkdir(gui.wd / gene, exist_ok=True)
        genemeta = data.metadata[gene]

        # do actual trimming
        for _id in data.gene_ids[gene]:
            record = genedata.pop(_id)
            record_meta = genemeta[_id]

            if _id in ignore_ids[gene]:
                record_meta['quality'] = 'manually dropped at trim_data'
                continue
            if not p.accept_rev and record_meta['is_rev']:
                record_meta['quality'] = 'disallowed reverse reads'
                continue
            try:
                record = trim_ends(record, p.min_phred, (p.trim_out, p.trim_of))
                record = mark_bad_stretches(record, p.min_phred, p.bad_stretch)
            except ValueError:
                record_meta['quality'] = 'low quality'
                continue
            except AttributeError:
                if not p.accept_nophred and 'accession' not in record_meta:
                    continue
                record_meta['quality'] = 'no phreds'
            # put back only the good records
            genedata[_id] = record

        # update dict of legal ids
        data.gene_ids[gene] = set(genedata.keys())

    # re-filter for shared entries
    shared_ids = set.intersection(*data.gene_ids.values())
    # re-loop seqdata
    LOG.debug('writing collated .fasta files')
    for gene, genedata in data.seqdata.items():
        # write to file
        with open(str(gui.wd / gene / (gene + '_all.fasta')), 'w') as fasta:
            SeqIO.write(genedata.values(), fasta, 'fasta')
        # write only records shared across all genes to fasta for MSA
        with open(str(gui.wd / gene / (gene + '.fasta')), 'w') as fasta:
            SeqIO.write([record for _id, record in genedata.items()
                         if _id in shared_ids], fasta, 'fasta')

        genemeta = data.metadata[gene]
        for _id in {_id for _id in genedata.keys() if _id not in shared_ids}:
            genemeta[_id]['quality'] = 'not in all genes'

    LOG.debug('writing metadata, deleting seqdata')
    shared.write_metadata(gui)
    data.seqdata.clear()

    shared.set_changed(gui, PAGE, False)
    if run_after:
        [run(gui) for run in run_after]
    return
