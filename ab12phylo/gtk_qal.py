# 2021 Leo Kaindl

import logging
import threading
from pathlib import Path
from time import sleep

import gi
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo import repo
from ab12phylo_cmd.filter import trim_ends, mark_bad_stretches
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
plt.set_loglevel('warning')
PAGE = 2


class qal_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.accept_rev.set_active(data.qal.accept_rev)
        iface.rev_handler = iface.accept_rev.connect('clicked', self.parse, None)
        iface.accept_nophred.set_active(True)
        iface.accept_nophred.connect('clicked', self.parse, None)

        iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=61, lower=0,
                                                      step_increment=1, page_increment=1))
        iface.min_phred.set_numeric(True)
        iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)

        for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
            wi = iface.__getattribute__(w_name)
            data.qal.__setattr__(w_name, int(wi.get_text()))
            wi.connect('changed', self.edit_numerical_entry)
            wi.connect('focus_out_event', self.parse)
        iface.min_phred.disconnect_by_func(self.parse)  # MARK leave this line alone

        for w_name in ['accept_rev', 'accept_nophred']:
            data.qal.__setattr__(w_name, iface.__getattribute__(w_name).get_active())

        for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
            wi.connect('key-press-event', self.edit_numerical_entry_up_down)

        # init row annotation
        iface.view_qal.set_model(data.qal_model)
        iface.view_qal.set_headers_visible(False)
        iface.view_qal.append_column(Gtk.TreeViewColumn(
            title='id', cell_renderer=Gtk.CellRendererText(),
            text=0, underline=2, strikethrough=3))
        iface.view_qal.set_tooltip_column(1)
        # crt = Gtk.CellRendererToggle(radio=False)
        # crt.props.indicator_size = 13
        # iface.view_qal.append_column(Gtk.TreeViewColumn(
        #     title='no phreds', active=1, cell_renderer=crt))

        sel = iface.view_qal.get_selection()
        sel.set_mode(Gtk.SelectionMode.MULTIPLE)
        sel.connect('changed', self.keep_visible, iface.view_qal,
                    iface.parallel_qal.props.vadjustment.props, iface.tempspace)
        iface.view_qal.connect(
            'check-resize', self.get_height_resize, iface.qal_spacer, [iface.qal_win])
        # in-preview deletion
        iface.view_qal.connect(
            'key_press_event', self.delete_and_ignore_rows, PAGE, sel)
        iface.qal_eventbox.connect_after(
            'button_press_event', self.select_seqs, PAGE, iface.zoomer,
            iface.view_qal, iface.tempspace)  # in-preview selection
        iface.qal_win.connect('scroll-event', self.xy_scale, PAGE)  # zooming

        # connect gene switcher
        iface.gene_handler = iface.gene_roll.connect('changed', self.select_gene_and_redo)

    def init_gene_roll(self):
        """
        Initialize gene switcher combo box with the
        previously selected gene from the dataset.
        """
        data = self.data
        iface = self.iface
        with iface.gene_roll.handler_block(iface.gene_handler):
            iface.gene_roll.remove_all()
            genes = list(data.genes)
            if not genes:
                return
            [iface.gene_roll.append_text(gene) for gene in genes]
            if len(genes) > 1:
                iface.gene_roll.insert_text(0, 'all')
                genes.insert(0, 'all')
            if data.gene_for_preview and data.gene_for_preview in genes:
                idx = genes.index(data.gene_for_preview)
            else:
                idx = 0
                data.gene_for_preview = 'all'
            iface.gene_roll.set_active(idx)

    def select_gene_and_redo(self, *args):
        """
        A page-independent handler for selecting a different gene.
        Currently, iface.gene_roll is only visible from one page, so a bit useless.
        :param args:
        :return:
        """
        page = self.iface.notebook.get_current_page()
        if page == 2:  # trim preview
            self.parse(self.iface.gene_roll, None)
            self.start_trim()

    def refresh(self):
        data = self.data
        iface = self.iface
        with iface.accept_rev.handler_block(iface.rev_handler):
            iface.accept_rev.set_active(data.search_rev and iface.accept_rev.get_active())
            iface.accept_rev.set_sensitive(data.search_rev)

        if not (self.wd / repo.PATHS.preview).exists() or 0 in data.qal_shape:
            self.start_trim()
            return

        # place the png preview
        self.load_image(iface.zoomer, PAGE, iface.qal_eventbox, self.wd / repo.PATHS.preview,
                        data.qal_shape[0] * self.get_hadj(), data.qal_shape[1])
        self.load_colorbar(iface.palplot)
        self.win.show_all()

    def reload_ui_state(self):
        for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
            self.iface.__getattribute__(w_name).set_text(str(self.data.qal.__getattribute__(w_name)))
        for w_name in ['accept_rev', 'accept_nophred']:
            self.iface.__getattribute__(w_name).set_active(self.data.qal.__getattribute__(w_name))

    def parse(self, widget, event):
        """
        Parse the content of a widget. Hitting enter in an entry still calls this, which is good.
        :param widget: The element to parse and inspect for changes
        :param event: passed by some signals -> focus_out_event
        :return:
        """
        data = self.data
        iface = self.iface
        LOG.debug('parsing for re-draw')
        w_name = widget.get_name()
        pre = None
        try:
            pre = data.qal.__getattribute__(w_name)
        except (KeyError, AttributeError):
            LOG.error('%s not in params' % w_name)
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
        self.delete_event(widget, event)

        # cause re-drawing if something changed
        if pre == now:
            LOG.debug('no change at parsing')
            return
        data.qal.__setattr__(w_name, now)
        if data.qal.trim_out > data.qal.trim_of:
            self.show_notification('cannot draw %d from %d' %
                                   (data.qal.trim_out, data.qal.trim_of))
            widget.set_text('0')
        else:
            self.set_changed(PAGE)

    def start_trim(self):
        """
        For non-empty gene set, start a background re-drawing thread
        and return to the main loop.
        :return:
        """
        data = self.data
        iface = self.iface
        # called after files are read
        if not data.genes or iface.thread.is_alive() \
                or iface.notebook.get_current_page() != PAGE:
            self.show_notification('Busy', secs=1)
            return

        if not data.seqdata:
            self.start_read(run_after=[self.start_trim])
            return

        data.gene_ids = {g: set(gd.keys()) for g, gd in data.seqdata.items() if g in data.genes}
        data.qal.accept_rev = iface.accept_rev.get_active()
        self.parse(iface.min_phred, None)  # annoying SpinButton

        LOG.debug('start-up redraw')
        data.qal_model.clear()
        sleep(.1)
        iface.thread = threading.Thread(target=self.do_trim)
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        # return to main loop

    def do_trim(self):
        """
        Iterate over records and trim, create a matrix representation
        of the valid characters and plot it.
        :return:
        """
        data = self.data
        iface = self.iface
        # parameters are up-to-date
        LOG.debug('re-draw trim preview')
        p = data.qal
        rows = list()
        data.gene_for_preview = p.gene_roll
        genes_for_preview = data.genes if p.gene_roll == 'all' else {p.gene_roll}
        ignore_ids = data.qal.ignore_ids if 'ignore_ids' in data.qal \
            else {g: set() for g in data.genes}  # regulated at delete_and_ignore_rows
        iface.text = 'creating matrix'
        LOG.debug(iface.text)
        iface.i = 0
        iface.k = sum([len(data.seqdata[gene]) for gene in genes_for_preview
                       if gene != 'no match']) + 4  # number of all records + extra
        shared_ids = set.intersection(*data.gene_ids.values())

        if not shared_ids:
            msg = 'No shared sequences between all genes (%s)' % ','.join(data.genes)
            self.show_notification(msg)
            LOG.warning(msg)
            sleep(.1)
            GObject.idle_add(self.stop_trim, True)
            return True

        try:
            for rid, gene in data.record_order:
                # skip records from other genes for the trimming preview
                if gene not in genes_for_preview or \
                        gene in ignore_ids and rid in ignore_ids[gene]:
                    continue

                iface.i += 1
                # maybe skip reversed seqs
                if not p.accept_rev and data.metadata[gene][rid]['is_rev']:
                    continue

                record = data.seqdata[gene][rid]
                try:
                    record = trim_ends(record, p.min_phred, (p.trim_out, p.trim_of), trim_preview=True)
                    record = mark_bad_stretches(record, p.min_phred, p.bad_stretch)
                    has_qal, is_bad = True, False
                    row = repo.seqtoint(record)
                except AttributeError:
                    # accept references anyway, but maybe skip no-phred ones
                    is_ref = 'accession' in data.metadata[gene][rid]
                    if not is_ref and not p.accept_nophred:
                        continue
                    has_qal, is_bad = is_ref, False
                    row = repo.seqtoint(record)
                except ValueError:
                    has_qal, is_bad = True, True
                    row = repo.seqtogray(record)
                rows.append(row)
                underline = not has_qal if rid in shared_ids else 4
                data.qal_model.append([rid, gene, underline, is_bad])

        except KeyError as ke:
            LOG.exception(ke)

        if not rows:
            msg = 'No sequence data remains'
            self.show_notification(msg)
            LOG.warning(msg)
            sleep(.1)
            GObject.idle_add(self.stop_trim)
            return True

        iface.text = 'tabularize'
        LOG.debug(iface.text)
        max_len = max(map(len, rows))
        array = np.array([row + repo.seqtoint(' ') * (max_len - len(row)) for row in rows])
        # make gaps transparent
        array = np.ma.masked_where(array > repo.toint('else'), array)
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
        f.set_figheight(array.shape[0] / repo.DPI)
        f.set_figwidth(array.shape[1] / repo.DPI)
        ax = f.add_subplot(111)
        ax.axis('off')
        f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        mat = ax.matshow(array, alpha=1, cmap=ListedColormap(repo.colors),
                         vmin=-.5, vmax=len(repo.colors) - .5, aspect='auto')

        iface.i += 1
        iface.text = 'save PNG'
        LOG.debug(iface.text)
        Path.mkdir(self.wd / repo.PATHS.preview.parent, exist_ok=True)
        f.savefig(self.wd / repo.PATHS.preview, transparent=True,
                  dpi=scale * repo.DPI, bbox_inches='tight', pad_inches=0.00001)
        plt.close(f)

        data.qal_shape[0] = array.shape[1]
        data.qal_shape[1] = self.get_height_resize(
            iface.view_qal, None, iface.qal_spacer, [iface.qal_win])

        if iface.rasterize.props.active:
            iface.text = 'place PNG'
            LOG.debug(iface.text)
            sleep(.05)
            self.load_image(iface.zoomer, PAGE, iface.qal_eventbox, self.wd / repo.PATHS.preview,
                            data.qal_shape[0] * self.get_hadj(), data.qal_shape[1])
        else:
            iface.text = 'place vector'
            LOG.debug(iface.text)
            canvas = FigureCanvas(f)  # a Gtk.DrawingArea
            canvas.set_size_request(data.qal_shape[0] * self.get_hadj(), data.qal_shape[1])
            try:
                ch = iface.qal_eventbox.get_child()
                if ch:
                    iface.qal_eventbox.remove(ch)
                iface.qal_eventbox.add(canvas)
            except Exception as ex:
                LOG.error(ex)
        iface.i += 1
        self.load_colorbar(iface.palplot)

        iface.text = 'idle'
        iface.frac = 1
        sleep(.1)
        GObject.idle_add(self.stop_trim)
        return True

    def stop_trim(self, errors=False):
        self.iface.thread.join()
        self.win.show_all()
        self.set_errors(PAGE, errors)
        self.iface.prog_bar.set_text('idle')
        LOG.info('qal thread idle')
        return False

    def trim_all(self, run_after=None):
        """
        Trim all SeqRecords in project_dataset.seqdata to sequence strings
        and write two collated .fasta files per gene, one only containing
        the records shared across all genes. If necessary, trace files are
        re-read before and will be deleted from memory afterwards for a
        smaller project file. Plot and place png previews of trimming result.
        :param run_after: [functions(gui)] when finished, usually flip page.
        :return:
        """
        data = self.data
        iface = self.iface

        if not data.seqdata:
            LOG.debug('re-reading files')
            self.start_read(run_after=[self.trim_all])
            return

        p = data.qal
        ignore_ids = p.ignore_ids if 'ignore_ids' in p \
            else {g: set() for g in data.genes}  # regulated at delete_and_ignore_rows
        LOG.debug('trim and filter all sequences')
        for gene, genedata in data.seqdata.items():
            Path.mkdir(self.wd / gene, exist_ok=True)
            genemeta = data.metadata[gene]

            # do actual trimming
            for _id in data.gene_ids[gene]:
                record = genedata.pop(_id)
                record_meta = genemeta[_id]

                if gene in ignore_ids and _id in ignore_ids[gene]:
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
                    is_ref = 'accession' in record_meta
                    if not p.accept_nophred and not is_ref:
                        continue
                    if not is_ref:
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
            with open(str(self.wd / gene / (gene + '_all.fasta')), 'w') as fasta:
                SeqIO.write(genedata.values(), fasta, 'fasta')
            # write only records shared across all genes to fasta for MSA
            with open(str(self.wd / gene / (gene + '.fasta')), 'w') as fasta:
                SeqIO.write([record for _id, record in genedata.items()
                             if _id in shared_ids], fasta, 'fasta')

            genemeta = data.metadata[gene]
            for _id in {_id for _id in genedata.keys() if _id not in shared_ids}:
                if 'quality' not in genemeta[_id]:
                    genemeta[_id]['quality'] = 'not in all genes'

        LOG.debug('writing metadata, deleting seqdata')
        self.write_metadata()
        data.seqdata.clear()

        self.set_changed(PAGE, False)
        self.set_changed(PAGE + 1, True)
        if run_after:
            [run() for run in run_after]
        return

    def delete_event(self, widget, event):
        return False
