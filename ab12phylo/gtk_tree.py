# 2021 Leo Kaindl

import itertools
import logging
import pickle
import threading
import xml.etree.ElementTree as ET
from copy import deepcopy
from math import ceil
from pathlib import Path
from time import sleep

import gi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import svgutils.transform as sg
import toyplot
import toytree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.cm import get_cmap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from toyplot import html, pdf, png, svg, color
from toytree.utils import TreeError, ToytreeError

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, Gdk, GLib

from ab12phylo import repo
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 7
UI_INFO = """
<ui>
  <popup name='context_menu'>
    <menuitem action='root' />
    <menuitem action='drop' />
    <menuitem action='collapse' />
    <menuitem action='select' />
    <menuitem action='extract' />
    <menuitem action='calculate' />
    <menuitem action='haplo_collapse' />
  </popup>
</ui>
"""


class tree_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface
        phy = data.phy

        # context menu
        action_group = Gtk.ActionGroup(name='context_menu')
        action_group.add_actions(
            [('root', None, 'Root here', None,
              'Root the tree at the selected node(s).', self.tree_modify),
             ('drop', None, 'Drop taxa', None,
              'Drop the selected nodes from the tree.', self.tree_modify),
             ('collapse', None, 'Collapse clade', None,
              'Collapse selected taxa to a leaf.', self.tree_modify),
             ('select', None, 'Select clade', None,
              'Select all descendants of the MRCA.', self.expand_popgen),
             ('extract', None, 'Extract clade', None,
              'Drop all other nodes from the tree.', self.tree_modify),
             ('calculate', None, 'Calculate diversity', None,
              'Calculate diversity metrics for the selected population.', self.start_popgen),
             ('haplo_collapse', None, 'Collapse haplotypes', None,
              'Collapse identical haplotyes in the selected population '
              'to a single leaf.', self.collapse_haplotypes)])

        uimanager = Gtk.UIManager()
        uimanager.add_ui_from_string(UI_INFO)
        self.win.add_accel_group(uimanager.get_accel_group())
        uimanager.insert_action_group(action_group)
        iface.contextmenu = uimanager.get_widget('/context_menu')

        for w_name in ['tree_eventbox', 'view_msa_ids', 'msa_eventbox']:
            iface.__getattribute__(w_name).connect(
                'button-press-event', self._on_right_click, iface.contextmenu)

        iface.query.connect('focus_out_event', self.select_popgen)
        iface.exclude.connect('focus_out_event', self.constrain_selection)

        iface.save_plot.connect('clicked', self.start_save)

        iface.view_msa_ids.set_model(data.tree_anno_model)
        crt = Gtk.CellRendererText(editable=True)
        crt.connect('edited', self._rename_taxon)
        col = Gtk.TreeViewColumn(title='id', markup=0, foreground=2, background=3,
                                 cell_renderer=crt)
        col.props.max_width = 80
        iface.view_msa_ids.append_column(col)
        iface.view_msa_ids.set_tooltip_column(1)

        iface.pop_metrics.set_model(data.pop_model)

        # set some plot_menu MenuButton images
        iface.flip.set_image(iface.flipim)
        iface.dist.set_image(iface.distim)
        iface.anno.set_image(iface.annoim)
        iface.filetypes.set_image(iface.saveim)
        iface.space.set_image(iface.spaceim)

        # connect zooming and selecting
        iface.tree_sel = iface.view_msa_ids.get_selection()
        iface.tree_sel.set_mode(Gtk.SelectionMode.MULTIPLE)
        iface.tree_sel.connect('changed', self.keep_visible, iface.view_msa_ids,
                               iface.parallel_tree.props.vadjustment.props,
                               iface.tempspace)
        iface.msa_eventbox.connect_after(
            'button_press_event', self.select_seqs, PAGE, iface.zoomer,
            iface.view_msa_ids, iface.tempspace)  # in-preview selection
        iface.tree_eventbox.connect_after(
            'button_press_event', self.select_seqs, PAGE, iface.zoomer,
            iface.view_msa_ids, iface.tempspace)  # in-preview selection
        iface.tree_pane.connect(
            'scroll-event', self.xy_scale, PAGE)  # zooming

        iface.tree_pane.set_position(599)
        self.iface.ratio = 600 / self.win.get_default_size()[0]
        self.win.connect_after('window-state-event', self._size_pane)

        iface.pane_handler = iface.tree_pane.connect_after(
            'notify::position', self._size_scrollbars,
            iface.tree_left_scrollbar, iface.tree_right_scrollbar,
            iface.tree_spacer, iface.view_msa_ids)

        iface.subtree.connect('clicked', self.expand_popgen)
        iface.calculate.connect('clicked', self.start_popgen)

        phy.tx = 'circular' if phy.circ else 'rectangular'
        phy.tx += '_TBE' if phy.tbe else '_FBP'
        phy.array = None
        iface.tree = None

    def reset_tree(self, *args):
        """
        Re-set the tree modifications by deleting
        the dictionary and the modified tree
        """
        if self.iface.thread.is_alive():
            self.iface.thread.join()
        self.data.phy.modify = dict()
        if 'g_lens' in self.data.phy:
            del self.data.phy.g_lens
        if (self.wd / repo.PATHS.modified_tree).is_file():
            (self.wd / repo.PATHS.modified_tree).unlink()
        self.start_phy()

    def refresh(self):
        """Re-view the page"""
        data = self.data
        iface = self.iface
        phy = data.phy
        self.reload_ui_state(page=PAGE)
        iface.tree_pane.set_position(600)

        if 'tx' in phy and (self.wd / phy.tx).with_suffix('.png').is_file() \
                and (self.wd / repo.PATHS.phylo_msa).is_file() \
                and 'shape' in phy:
            try:
                # load tree PNG
                self.load_image(iface.zoomer, PAGE, iface.tree_eventbox,
                                self.wd / (phy.tx + '.png'), h=phy.shape[1])
                # load MSA PNG
                self.load_image(iface.zoomer, PAGE, iface.msa_eventbox,
                                self.wd / repo.PATHS.phylo_msa,
                                w=phy.shape[0] * self.get_hadj() * 2, h=phy.shape[1])
                # load colorbars
                self.load_colorbar(iface.palplot1)
                if phy.axis and len(data.genes) > 1:
                    self.load_colorbar(iface.gbar, gbar=True)
                    iface.gbar.set_visible(True)
                self.win.show_all()
                return
            except GLib.Error:
                pass

        self.start_phy()
        self.win.show_all()

    def select_popgen(self, widget, ev):
        """Select samples for popgen calculations by name / species matching"""
        data = self.data
        iface = self.iface
        phy = data.phy
        phy.query = widget.get_text()
        phy.exclude = iface.exclude.get_text()

        if not phy.query and not phy.exclude:
            return
        mo = data.tree_anno_model

        search_space = [_id + sp for _id, sp in mo.get_column((0, 1))]
        sel = list()

        # first match
        if not phy.query:
            sel = list(range(len(mo)))
        else:
            sel = [i for i, label in enumerate(search_space)
                   if repo.rgx(phy.query).match(label)]

        # then exclude
        if phy.exclude:
            sel = [i for i, label in enumerate(search_space) if
                   i in sel and not repo.rgx(phy.exclude).match(label)]

        # then select in the model
        iface.tree_sel.unselect_all()
        [iface.tree_sel.select_path(i) for i in sel]

    def constrain_selection(self, widget, ev, passed=None):
        data = self.data
        iface = self.iface
        phy = data.phy
        if passed:
            sel_idx = passed
        else:
            sel_idx = [tp.get_indices()[0] for tp
                       in iface.tree_sel.get_selected_rows()[1]]
        if not sel_idx:
            return

        # get exclusion search space
        idx_and_label = [(idx, _id + sp) for idx, (_id, sp) in enumerate(
            data.tree_anno_model.get_column((0, 1))) if idx in sel_idx]

        # then exclude
        phy.exclude = iface.exclude.get_text()
        sel_idx = [i for i, label in idx_and_label if
                   not repo.rgx(phy.exclude).match(label)]
        if widget:
            # then select in the model
            iface.tree_sel.unselect_all()
            [iface.tree_sel.select_path(i) for i in sel_idx]
        else:
            return sel_idx

    def expand_popgen(self, widget, sel_idx=None):
        """Select a subtree starting from a GtkTreeSelection"""
        data = self.data
        iface = self.iface
        phy = data.phy
        if not sel_idx:
            sel_idx = [tp.get_indices()[0] for tp
                       in iface.tree_sel.get_selected_rows()[1]]
        if not sel_idx:
            self.show_notification('Cannot expand to subtree, no selection', secs=1)
            return
        if not self.iface.tree:
            self.show_notification('No tree in memory, re-drawing first', secs=5)
            self.start_phy(run_after=(self.expand_popgen, (None, sel_idx)))
            return

        in_blue = lambda c: c.replace('style="normal"', 'style="fill:%s"'
                                      % repo.tohex(repo.rocket[0]))

        # get a list of all tip labels and whether they are selected
        id_and_ifsel = [(in_blue(_id), i in sel_idx) for i, _id in enumerate(
            iface.tree.get_tip_labels()[::-1])]

        # move up in the tree: find the MRCA
        mrca = iface.tree.get_mrca_idx_from_tip_labels(
            names=[_id for _id, sel in id_and_ifsel if sel])

        # move down in the tree: get tip positions
        sel_idx = [[_id for _id, sel in id_and_ifsel].index(_id)
                   for _id in iface.tree.get_tip_labels(mrca)]

        if phy.exclude:
            sel_idx = self.constrain_selection(None, None, sel_idx)

        [iface.tree_sel.select_path(i) for i in sel_idx]
        iface.pop_size.set_text('population size: %d' % len(sel_idx))
        iface.pop_size.set_visible(True)

    def start_popgen(self, widget=None, haplo_collapse=False):
        data = self.data
        iface = self.iface
        phy = data.phy
        data.pop_model.clear()
        self.init_popgen_columns(iface.pop_metrics)
        sleep(.1)  # prevent crash on re-load->Ctrl+A->calculate

        mo, tps = iface.tree_sel.get_selected_rows()
        tps = [tp[0] for tp in tps]
        if haplo_collapse:
            tps = haplo_collapse

        if len(tps) < 2:
            self.show_notification('invalid selection', secs=1)
            return

        if phy.array is None:
            self.start_phy()
            return

        for w_name in ['gap_share', 'unk_share']:
            phy.__setattr__(w_name, iface.__getattribute__(
                w_name).get_adjustment().props.value)

        # re-direct to thread
        iface.i = 0
        iface.k = len(data.genes) + 2
        iface.text = 'calculating'
        iface.thread = threading.Thread(target=self.do_popgen, args=[tps, haplo_collapse])
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()

    def do_popgen(self, tps, haplo_collapse):
        data = self.data
        iface = self.iface
        phy = data.phy
        iface.pop_size.set_text('population size: %d' % len(tps))
        iface.pop_size.set_visible(True)

        # show the ids of the selected samples
        buf = Gtk.TextBuffer()
        buf.props.text = ', '.join(
            [_id for i, _id in enumerate(self.data.tree_anno_model
                                         .get_column(0)) if i in tps])
        iface.pop_ids.set_buffer(buf)

        # crop the array and shift the entries
        ar = np.flipud(phy.array)[tps, :]
        # shift 4|>7 to 5 and 6|7 to 4
        ar = np.array([[5 if x == 4 or x > 7 else (4 if x in {6, 7} else x)
                        for x in _ar] for _ar in ar])
        start, overall, drops, repeats = 0, set(), set(), dict()
        for i, gene in enumerate(data.genes):
            iface.text = 'popgen: %s' % gene
            LOG.debug(iface.text)
            _range = range(start, start + phy.g_lens[gene])
            overall.update(_range)
            row, drops, repeats = _per_gene_diversity(gene, phy, ar, _range)
            data.pop_model.append(row)
            start += phy.g_lens[gene] + len(repo.SEP)
            iface.i += 1

        if len(data.genes) > 1:
            iface.text = 'popgen: overall'
            LOG.debug(iface.text)
            row, drops, repeats = _per_gene_diversity('overall', phy, ar, overall)
            data.pop_model.append(row)
            iface.i += 1

        # show the table
        iface.tree_stack.set_visible_child_name('popgen')

        self._highlight_duplicates(tps, repeats)

        if haplo_collapse:
            self.collapse_haplotypes(None, return_loop=True)
        else:
            # plot using tps and drops
            self._do_phy5((tps, drops, repeats))

    def _highlight_duplicates(self, tps, repeats):
        """
        Highlight duplicate sequences by shading the ID column in matching color.
        :param tps: indices [int] of selected tree tips
        :param repeats: {idx: [occurrences of seq found at idx, given it occurs several times]}
        :return:
        """
        mo = self.data.tree_anno_model
        # remove old color markup
        for row in mo:
            if row[3] != self.iface.AQUA:
                row[3] = None

        colors = [repo.tohex(c[:3]) for c in get_cmap('Set1', len(repeats)).colors]
        for color, idcs in zip(colors, repeats.values()):
            for idx in idcs:
                mo[tps[idx]][3] = color

    def collapse_haplotypes(self, action, return_loop=False, sel_idx=None):
        if not return_loop:
            if not sel_idx:
                sel_idx = [tp.get_indices()[0] for tp
                           in self.iface.tree_sel.get_selected_rows()[1]]
            if not sel_idx:
                self.show_notification('Cannot %s, no selection' % action.get_name(), secs=2)
                return
            if not self.iface.tree:
                self.show_notification('No tree in memory, re-drawing first', secs=5)
                self.start_phy(run_after=(self.collapse_haplotypes, (action, False, sel_idx)))
                return
            # first re-find the haplotyes
            self.start_popgen(haplo_collapse=sel_idx)
        else:
            # then collapse
            props = self.iface.spacespin.get_adjustment().props
            props.value = max(2.4, props.value)
            mod = self.data.phy.modify
            candidates = dict()
            for i, (n, sp, c) in enumerate(self.data.tree_anno_model.get_column((0, 1, 3))):
                if not c or c == self.iface.AQUA:
                    continue
                if sp.startswith('REF_'):
                    n = sp.split(':')[0]
                if c in candidates:
                    candidates[c][i] = n  # tips[i]
                else:
                    candidates[c] = {i: n}  # tips[i]}
            if candidates:
                if 'rename' not in mod:
                    mod['rename'] = dict()
                if 'drop' not in mod:
                    mod['drop'] = list()
            for haplo_dict in candidates.values():
                # rename the remaining node
                mod['rename'][haplo_dict[min(haplo_dict)]] = '~'.join(
                    sorted([t for t in haplo_dict.values() if not t.startswith('REF_')]) + \
                    sorted([t for t in haplo_dict.values() if t.startswith('REF_')]))
                # delete the others
                mod['drop'].extend([v for k, v in haplo_dict.items() if k != min(haplo_dict)])

            self.iface.run_after = self.start_phy, ()

            sleep(.1)
            GObject.idle_add(self._stop_phy)
            return True

    def _rename_taxon(self, cell, path, new_text):
        mo = self.data.tree_anno_model
        old_text = mo[path][0]
        if mo[path][1].startswith('REF_'):
            old_text = mo[path][1].split(':')[0]
        if old_text == new_text:
            return
        elif new_text in mo.get_column(0) + [i.split(':')[0] for i in mo.get_column(1)]:
            return
        LOG.debug('renaming %s to %s' % (old_text, new_text))
        mo[path][0] = new_text
        new_text = new_text.replace(', ', '~').replace(',', '~')
        old_text = old_text.replace(', ', '~').replace(',', '~')
        mod = self.data.phy.modify
        if 'rename' not in mod:
            mod['rename'] = {old_text: new_text}
        else:
            # save twice
            for k, v in mod['rename'].items():
                if v == old_text:
                    mod['rename'][k] = new_text
            mod['rename'][old_text] = new_text
        LOG.debug(mod['rename'])

    def init_popgen_columns(self, tv):
        if tv.get_n_columns() > 0:
            return
        for i, ti in enumerate(['gene', 'sites', 'S', 'k', 'π', 'Watterson\'s θ',
                                'Tajima\'s D', 'unique seqs', 'gaps', 'unknown']):
            col = Gtk.TreeViewColumn(
                title=ti, cell_renderer=Gtk.CellRendererText(), text=i)
            col.set_resizable(True)
            tv.append_column(col)
        tv.realize()
        tv.columns_autosize()

    def init_sp_genes(self, sp_genes, genes, sel_gene=None):
        """
        Init the selector for which gene to draw species annotations from
        :param sp_genes: the GtkComboBoxText to display the genes
        :param genes: the iterable of genes in the analysis
        :param sel_gene: an optional selected gene
        :return: sel_gene: if not set before, is 'best pid' now
        """
        sp_genes.remove_all()
        genes = list(genes)
        if len(genes) > 1:
            genes.insert(0, 'best pid')
        [sp_genes.append_text(gene) for gene in genes]
        if sel_gene and sel_gene in genes:
            idx = genes.index(sel_gene)
        else:
            idx = 0
            sel_gene = 'best pid'
        sp_genes.set_active(idx)
        return sel_gene

    def start_save(self, *args):
        """Handle the "Export the current plot now" save functionality via a thread."""
        iface = self.iface
        phy = self.data.phy
        if 'canvas' not in iface:
            self.start_phy()
            return

        # parse plot settings only when needed
        [phy.__setattr__(w_name, iface.__getattribute__(
            w_name).get_active()) for w_name in
         ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
          'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]
        phy.flip = iface.flipspin.props.adjustment.props.value
        phy.dist = iface.distspin.props.adjustment.props.value
        phy.space = iface.spacespin.props.adjustment.props.value
        phy.sel_gene = iface.sp_genes.get_active_text()

        # re-direct to thread
        iface.thread = threading.Thread(target=self._do_save)
        iface.k = 2 + 2 * phy.pdf + 2 * phy.svg + phy.pmsa + phy.html
        iface.i = 0
        iface.text = 'read tree'
        sleep(.1)
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        return
        # return to main loop

    def _do_save(self, *args):
        iface = self.iface
        self.save_tree()
        self.save_msa()
        self.iface.text = 'idle'
        sleep(.1)
        GObject.idle_add(self._stop_save)
        return True

    def _stop_save(self, *args):
        self.iface.thread.join()
        LOG.debug('phylo saving idle')
        self.win.show_all()
        self.data.change_indicator[PAGE] = False
        return False

    def save_tree(self, *args):
        """Save the current tree in different formats"""
        data = self.data
        iface = self.iface
        phy = data.phy
        LOG.debug('rendering tree')
        if 'canvas' not in iface:
            self.start_phy()
            return
        tx = self.wd / phy.tx
        iface.text = 'rendering png'
        png.render(iface.canvas, str(tx.with_suffix('.png')), scale=1.6)
        iface.i += 1
        if phy.pdf:
            iface.text = 'rendering pdf'
            pdf.render(iface.canvas, str(tx.with_suffix('.pdf')))
            iface.i += 1
        if phy.svg:
            iface.text = 'rendering svg'
            svg.render(iface.canvas, str(tx.with_suffix('.svg')))
            iface.i += 1
        if phy.html:
            iface.text = 'rendering html'
            html.render(iface.canvas, str(tx.with_suffix('.html')))
            iface.i += 1

    def save_msa(self, popgen=None):
        """Save the current msa in different formats"""
        data = self.data
        iface = self.iface
        phy = data.phy
        LOG.debug('rendering msa')
        Path.mkdir(self.wd / repo.PATHS.phylo_msa.parent, exist_ok=True)
        msa = self.wd / repo.PATHS.phylo_msa.parent / repo.PATHS.phylo_msa.stem
        dpi = iface.scale * repo.DPI
        iface.text = 'rendering png'
        iface.fig.savefig(msa.with_suffix('.png'), transparent=True, dpi=dpi)
        iface.i += 1

        if popgen:
            return
        if phy.svg:
            iface.text = 'rendering svg with msa -- SLOW'
            iface.fig.savefig(msa.with_suffix('.svg'), transparent=True)
            iface.i += 1
        if phy.pdf:
            iface.text = 'rendering pdf with msa -- SLOW'
            iface.fig.savefig(msa.with_suffix('.pdf'), transparent=True, dpi=dpi)
            iface.i += 1

        if phy.pmsa:
            if phy.svg:
                iface.text = 'collating svg tree and msa -- SLOW'
                LOG.info(iface.text)
                tree = sg.fromfile((self.wd / phy.tx).with_suffix('.svg'))
                m = sg.fromfile(msa.with_suffix('.svg'))
                mw, mh = [float(i[:-2]) * 96 / 72 for i in m.get_size()]  # now in px
                pt = tree.getroot()
                tw, th = [float(i[:-2]) for i in tree.get_size()]  # already in px
                pm = m.getroot()
                pm.moveto(tw, 0, th / mh / 6, th / mh * 96 / 72)
                fig = sg.SVGFigure()  # '%.0fpx' % (tw + th / mh / 8), '%.0fpx' % th)
                fig.append([pt, pm])
                fig.save(self.wd / (phy.tx + '_msa.svg'))
            else:
                iface.text = 'rendering png with msa'
                LOG.debug(iface.text)
                old_wi = iface.canvas.width
                wi = phy.array.shape[1] * 3
                im = Image.open(msa.with_suffix('.png'))
                im = im.convert('RGB')
                im = im.resize((wi, int(iface.canvas.height)), resample=Image.NEAREST)
                iface.canvas.width += wi
                im_axes = iface.canvas.image(im, bounds=(old_wi, old_wi + wi - 20,
                                                         0, iface.canvas.height))
                png.render(iface.canvas, str(self.wd / (phy.tx + '_msa.png')), scale=1.6)
            iface.i += 1

    def start_phy(self, modify=None, run_after=None):
        """Get settings and block GUI"""
        data = self.data
        iface = self.iface
        phy = data.phy
        if iface.thread.is_alive():
            self.show_notification('Busy', secs=1)
            return
        LOG.debug('start_plotting')
        if modify:
            k, v = modify.popitem()
            if k == 'root':  # replace old root
                phy.modify[k] = v
            else:
                phy.modify[k] = phy.modify.get(k, list()) + v

        # empty placeholder early to avoid over-high graphics
        for wi in [iface.tree_eventbox, iface.msa_eventbox]:
            ch = wi.get_child()
            if ch:
                ch.destroy()

        # parse plot settings only when needed
        [phy.__setattr__(w_name, iface.__getattribute__(
            w_name).get_active()) for w_name in
         ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
          'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]
        phy.flip = iface.flipspin.props.adjustment.props.value
        phy.dist = iface.distspin.props.adjustment.props.value
        phy.space = iface.spacespin.props.adjustment.props.value
        phy.sel_gene = iface.sp_genes.get_active_text()

        # re-direct to thread
        iface.thread = threading.Thread(target=self._do_phy1)
        iface.k = (phy.rect + phy.circ + phy.unro) * 2 + 9 \
                  + 2 + phy.html + phy.pmsa + 2 * phy.svg + 2 * phy.pdf  # saving
        iface.i = 0
        iface.text = 'read tree'
        sleep(.1)
        iface.run_after = run_after
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        return
        # return to main loop

    def _do_phy1(self):
        data = self.data
        iface = self.iface
        phy = data.phy
        drops = set()

        # pick a tree file
        if phy.modify and (self.wd / repo.PATHS.modified_tree).is_file():
            LOG.debug('using modified tree')
            tree_file = self.wd / repo.PATHS.modified_tree
        else:
            tree_file = self.wd / repo.PATHS.fbp if phy.fbp \
                else self.wd / repo.PATHS.tbe
            if not tree_file.is_file():
                GObject.idle_add(self._stop_phy, 'no tree found at %s' % str(tree_file))
                sleep(.1)
                return True

        # read tree file
        tr = toytree.tree(open(tree_file, 'r').read(), tree_format=0)

        if 'collapse' in phy.modify:
            names = [i for i in phy.modify['collapse'] if i in tr.get_tip_labels()]
            if names:
                # move references to the back
                names = sorted([t for t in names if not t.startswith('REF_')]) \
                        + sorted([t for t in names if t.startswith('REF_')])
                try:
                    mrca_idx = tr.get_mrca_idx_from_tip_labels(names=names)
                    mrca = tr.idx_dict[mrca_idx]
                    drops = set(mrca.get_leaf_names()) - set(names)
                    [t.detach() for t in mrca.get_children()]
                    mrca.name = '~'.join(names)
                    phy.modify['collapse'] = [' ']  # only remember that, not what
                    tr = toytree.tree(tr.newick)
                except (TreeError, ToytreeError) as te:
                    GObject.idle_add(self._do_phy2, [str(te)])
                    sleep(.1)
                    return True

        if 'drop' in phy.modify:
            names = [i for i in phy.modify['drop'] if i in tr.get_tip_labels()]
            if names:
                tr = tr.drop_tips(names=names)
                tr = toytree.tree(tr.newick)

        if 'root' in phy.modify:
            try:
                names = [i for i in phy.modify['root'] if i in tr.get_tip_labels()]
                if names:
                    tr = tr.root(names=names)
            except ToytreeError as te:
                GObject.idle_add(self._do_phy2, [str(te)])
                sleep(.1)
                return True

        if 'extract' in phy.modify:
            names = [i for i in phy.modify['extract'] if i in tr.get_tip_labels()]
            if names:
                treenode = tr.idx_dict[tr.get_mrca_idx_from_tip_labels(names=names)]
                tr = toytree.tree(treenode.write())

        if 'rename' in phy.modify:
            for old, new in phy.modify['rename'].items():
                try:
                    tr.idx_dict[tr.get_mrca_idx_from_tip_labels(names=[old])].name = new
                except (TreeError, ToytreeError) as te:
                    pass

        if phy.modify:
            # save the modified tree
            tr.write(self.wd / repo.PATHS.modified_tree, tree_format=0)

        # fetch tips labels, now top-to-bottom
        phy.tips = list(reversed(tr.get_tip_labels()))
        separated_tips = {item for t in phy.tips for item in t.split('~')}
        iface.i += 1
        iface.tree = tr

        iface.text = 'read metadata'
        # read in tabular data
        df = pd.read_csv(self.wd / repo.PATHS.tsv, sep='\t', dtype={'id': str})
        df.set_index('id', inplace=True)
        iface.tempspace.df = df
        df = df.copy(deep=True)

        r, b, e = 'reference_species', 'BLAST_species', 'extra_species'
        if r in df:
            # overwrite BLAST with reference
            df.loc[df[r].notna(), b] = df.loc[df[r].notna(), r]
        if b in df:
            df = df.rename(columns={b: 'species'})
            phy.did_BLAST = True
        # crop to relevant columns
        df = df[[c for c in ['gene', 'accession', 'species', 'pid'] if c in df]]

        # carry over renamed taxa to metadata entries
        if 'rename' in phy.modify:
            df = df.reset_index()
            for old_id, new_id in phy.modify['rename'].items():
                dd = df.loc[(df['id'] == old_id)].copy()
                dd.loc[:, 'id'] = new_id
                df = df.append(dd)
            df.set_index('id', inplace=True)

        df = df.fillna(value={'pid': 0})
        pid = 'pid' in df
        if phy.sel_gene not in [None, 'best pid'] or len(data.genes) == 1:
            # picked species labels for a single, unique gene
            df = df.loc[df['gene'] == phy.sel_gene]
            df = df.loc[[t for t in separated_tips if t in df.index]]  # skip dropped nodes
            phy.ndf = df.to_dict(orient='index')
        else:
            # picked species labels from best pid across several genes
            # split into per-gene dataframes, filtering out filtered out genes
            dfs = [i[1] for i in df.groupby(df['gene']) if i[0] in data.genes]
            # find row with highest pid from all genes
            phy.ndf = dict()
            for t in separated_tips:
                if any([t not in d.index for d in dfs]):
                    continue  # a dropped node
                rows = [d.loc[t] for d in dfs]
                if pid:
                    phy.ndf[t] = [dict(r) for r in rows if r.pid ==
                                  max([q.pid for q in rows])][0]
                else:
                    phy.ndf[t] = [dict(r) for r in rows][0]

        # fetch species labels for compound nodes
        for compound_tip in set(phy.tips) - separated_tips:
            df = pd.DataFrame.from_records(
                phy.ndf[t] for t in compound_tip.split('~') if t in phy.ndf)
            if 'accession' in df:
                # extract simple species label
                df.loc[df['accession'].notna(), 'species'] = \
                    df.loc[df['accession'].notna()].species.replace(' strain.+', '', regex=True)
            if pid:
                df = df.loc[df.pid == df.pid.max()]
                # split compound species, then re-collate
                species = ' / '.join({item for sp in set(df.species) for item in sp.split(' / ')})
                phy.ndf[compound_tip] = {'species': species, 'pid': df.pid.max()}
            else:
                phy.ndf[compound_tip] = dict()

        iface.i += 1
        GObject.idle_add(self._do_phy2, None, drops)
        sleep(.1)
        return True

    def _do_phy2(self, errors=None, drops=None):
        """Bounce to idle to size the id column"""
        data = self.data
        iface = self.iface
        phy = data.phy
        iface.thread.join()
        in_blue = '<span style="normal">%s</span>'

        if errors:
            iface.text = 'idle'
            self.show_notification('Tree modification failed', errors)
            return
        if drops:
            self.show_notification('Collapsing a branch removed '
                                   'taxa that were not selected', drops, secs=6)

        iface.text = 'fill id column'
        LOG.debug(iface.text)
        data.tree_anno_model.clear()

        # write samples to model
        for t in phy.tips:
            if t not in phy.ndf and '~' not in t:
                data.tree_anno_model.append([t, '', None, iface.AQUA])
                continue  # a dropped node
            nd = phy.ndf[t]
            # make sure pid field is present
            pid = nd.get('pid', 0)
            if pd.isna(pid):
                pid = 0
            nd['pid'] = pid
            # make sure species field is present
            sp = nd.get('species', '')
            if pd.isna(sp):
                sp = ''
            nd['species'] = sp
            if 'accession' in nd and not pd.isna(nd['accession']):
                sp = t + ': ' + sp
                if len(data.genes) > 1 and 'strain' in nd['species']:
                    t = nd['species'][nd['species'].find('strain') + 7:]
                else:
                    t = nd['accession']
                c = iface.BLUE
            else:
                c = None  # iface.FG

            # edit compound tips
            if '~' in t:
                t = t.split('~')
                # for i, tp in enumerate(t):
                #     if tp.startswith('REF_'):
                #         nd = phy.ndf[tp]
                #         if len(data.genes) > 1 and 'strain' in nd['species']:
                #             t[i] = in_blue % nd['species'].split(' strain ')[1]
                #         else:
                #             t[i] = in_blue % nd['accession']
                t = ', '.join(t)
            data.tree_anno_model.append([t, sp, c, None])

        iface.view_msa_ids.set_margin_bottom(50 * phy.axis)
        iface.view_msa_ids.columns_autosize()
        iface.view_msa_ids.realize()
        sleep(.1)

        # re-direct to thread
        iface.thread = threading.Thread(target=self._do_phy3)
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()

    def _do_phy3(self):
        """Write new, annotated MSA, annotated nwk tree and plot tree"""
        data = self.data
        iface = self.iface
        phy = data.phy
        iface.text = 'write new MSA'
        errors = list()

        # read MSA
        records = {r.id: r.seq for r in SeqIO.parse(self.wd / repo.PATHS.msa, 'fasta')}
        # carry over renamed taxa to seqdata
        if 'rename' in phy.modify:
            for old_id, new_id in phy.modify['rename'].items():
                if old_id in records:
                    records[new_id] = records[old_id]
        phy.seqdata = list()
        # write re-named re-ordered MSA
        with open(self.wd / repo.PATHS.msa_anno, 'w') as new_msa:
            for tip in phy.tips:
                if tip not in records and '~' not in tip:
                    phy.seqdata.append(None)
                    continue  # a dropped node
                elif '~' in tip:
                    try:
                        ar = np.array([repo.seqtoint(str(records[t]))
                                       for t in tip.split('~')])
                    except KeyError as ke:
                        LOG.exception(ke)
                        GObject.idle_add(self.reset_tree)
                        sleep(.1)
                        return True
                    seq = repo.inttoseq([np.argmax(np.bincount(ar[:, i])) for i in range(ar.shape[1])])
                    records[tip] = seq

                if repo.SEP in records[tip]:
                    # get seq and cut out artificial separator
                    seq = str(records[tip]).replace(repo.SEP, '')
                elif 'g_lens' in phy:
                    # introduce separators for the sake of plotting
                    seq = str(records[tip])
                    pos = 0
                    for i in range(len(data.genes) - 1):
                        nlen = phy.g_lens[data.genes[i]]
                        seq = seq[:pos + nlen] + repo.SEP + seq[pos + nlen:]
                        pos += nlen + len(repo.SEP)
                else:
                    seq = str(records[tip])

                # get metadata
                try:
                    entry = phy.ndf[tip]
                except KeyError as ke:
                    LOG.exception(ke)
                    GObject.idle_add(self.reset_tree)
                    sleep(.1)
                    return True

                des = entry.get('species', '')
                pid = ('%.2f' % entry['pid']).rstrip('0').rstrip('.')
                if des and pid != '0':
                    des += '; pid %s' % pid

                phy.seqdata.append(repo.seqtoint(records[tip]))

                # rename REFs
                if tip.startswith('REF_') and '~' not in tip:
                    if len(data.genes) > 1 and 'strain' in entry['species']:
                        # # <Strain> <Species> strain <Strain>; pid <pid>
                        # tip = entry['species'].split(' strain ')[1]
                        # <Species> strain <Strain>; pid <pid>
                        tip = des
                        des = ''
                    else:
                        tip = entry['accession']

                SeqIO.write(SeqRecord(Seq(seq), id=tip, description=des), new_msa, 'fasta')
        LOG.debug('wrote updated MSA')
        iface.i += 1

        # get lengths of genes
        if 'g_lens' not in phy:
            phy.g_lens = {gene: len(seq) for gene, seq in
                          zip(data.genes, records[next(iter(records))].split(repo.SEP))}
        multi = True if len(data.genes) > 1 else False

        # cat: REF_ -> 1, no_BLAST_hit -> -1, regular -> 0
        # pid: discrete index in rocket palette
        # confident: 1 if >= flip threshold -> blue, else 0 -> red

        iface.text = 'iterate tree: 1'
        cat_ignore = -1
        in_blue = '<span style="fill:%s">%s</span>' % (repo.tohex(repo.rocket[0]), '%s')

        # guess if support val.s are in [0-1] or in [0-100], then scale with 1 or 1/100 resp.
        sup_scaler = 1 if 1.1 > max([float(s) for s in
                                    iface.tree.get_node_values('support') if s]) else 100
        # rename reference nodes and add features
        for node in iface.tree.treenode.traverse():
            if not node.support:
                node.add_feature('support', 0)
            elif type(node.support) != float:
                node.support = float(node.support)
            if node.is_leaf() and node.name in phy.ndf:
                entry = phy.ndf[node.name]

                if node.name.startswith('REF_'):
                    if multi and 'strain' in entry['species']:
                        entry['species'], strain = entry['species'].split(' strain ')
                        phy.ndf['strain'] = phy.ndf[node.name]
                        node.name = strain
                    else:
                        if 'accession' in entry:
                            acc = entry['accession']
                            phy.ndf[acc] = phy.ndf[node.name]
                            node.name = acc
                        else:  # collapsed node of refs
                            names = node.name.split('~')
                            for i, t in enumerate(names):
                                if t.startswith('REF_'):
                                    ndf = phy.ndf[t]
                                    if multi and 'strain' in ndf['species']:
                                        names[i] = in_blue % ndf['species'].split(' strain ')[1]
                                    else:
                                        names[i] = phy.ndf[t]['accession']
                            node.name = ', '.join(names)
                    node.add_feature('cat', 1)  # blue
                    node.add_feature('species', entry['species'])
                    node.add_feature('pid', 0)  # this is the index of blue in
                    # the rocket palette, which is used to display BLAST pid

                elif 'species' in entry:
                    if '~' in node.name:
                        names = node.name.split('~')
                        for i, t in enumerate(names):
                            if t.startswith('REF_'):
                                ndf = phy.ndf[t]
                                if multi and 'strain' in ndf['species']:
                                    names[i] = in_blue % ndf['species'].split(' strain ')[1]
                                else:
                                    names[i] = in_blue % ndf['accession']
                        node.name = ', '.join(names)
                    node.add_feature('species', entry['species'])
                    node.add_feature('cat', 0 if entry['species'].strip()
                                                 or not phy.spec else 2)  # black else dark_red
                    node.add_feature('pid', 1 + ceil(min(100 - entry['pid'], 20) / 2))  # map to the palette

                else:
                    raise ValueError('No species in ndf dataframe should be impossible')
                    # node.add_feature('cat', 0 if not phy.did_BLAST else 2)  # black else dark_red
                    # node.add_feature('pid', 1 if not phy.did_BLAST else 11)
                    # dark gray if no BLAST, else bad light color signifying no hit

                # meaningless for tips, but to make the fields exist
                node.add_feature('size', 0)
                node.add_feature('confident', 0)
            else:
                # use support values for node size
                node.support /= sup_scaler
                node.add_feature('species', '')
                node.add_feature('confident', 1 if node.support >= phy.flip else 0)
                node.add_feature('pid', -1)  # dark_red for no BLAST hit
                node.add_feature('cat', cat_ignore)  # black, just so the field exists

        iface.i += 2

        if phy.nwk:
            iface.text = 'export nwk'
            exp_tree = deepcopy(iface.tree)
            for node in exp_tree.treenode.traverse():
                if node.is_leaf():
                    try:
                        if node.name in node.species:
                            node.name = node.species
                        else:
                            node.name += ' ' + node.species
                            node.name += ' %.2f' % node.pid
                    except AttributeError:
                        continue
            # write new-ish newick file
            tf = self.wd / repo.PATHS.fbpn if phy.fbp else self.wd / repo.PATHS.tben
            exp_tree.write(tf, tree_format=0)

        # tree plotting tweaks
        iface.text = 'iterate tree: 2'
        for node in iface.tree.treenode.traverse():
            node.dist = max(node.dist, phy.dist)
            if phy.supp:
                node.support = int(round(node.support * 100))
        iface.i += 1

        iface.text = 'prepare plot attributes'
        iface.view_msa_ids.set_margin_bottom(50 * phy.axis)
        iface.view_msa_ids.realize()
        sleep(.1)
        # get model height
        phy.height = iface.view_msa_ids.get_allocated_height() + 50 * phy.axis

        # pre-calculate plot attributes
        # internal node size -> support
        min_size = .2
        ns = [min_size + s / 7.5 for s in iface.tree.get_node_values('support', 1, 1)]

        # internal node color -> support
        red_blue = [color.rgb(*c) for c in [repo.red, repo.blue]]
        nc = [red_blue[i] for i in iface.tree.get_node_values('confident', 1, 1)]

        # primary tip label color -> category
        cat_rgb = [color.rgb(*c) for c in [repo.black, repo.blue, repo.dark_red]]
        tlc = [cat_rgb[n] for n in list(iface.tree.get_node_values(
            'cat', 1, 1)) if n != cat_ignore][::-1]

        # secondary tip label color -> BLAST pid on species
        tl = iface.tree.get_tip_labels()
        if phy.spec:
            # rocket_rgb = [color.rgb(*c) for c in rocket]
            rocket_css = [repo.tohex(c) for c in repo.rocket]
            tl = ['%s <i><span style="fill:%s">%s</span></i>'
                  % (t.name, rocket_css[t.pid], t.species) for t in
                  iface.tree.treenode.traverse('postorder') if t.is_leaf()][::-1]

        # plot tree and msa
        if phy.circ:
            phy.tx = 'circular'
            ly = 'c'
            w, h = 1600, 1600
            iface.canvas = toyplot.Canvas(width=w, height=h)
            iface.axes = iface.canvas.cartesian(bounds=(200, w - 200, 200, h - 200))
        elif phy.rect:
            phy.tx = 'rectangular'
            ly = 'r'
            # get dim for canvas
            w, h = 1200 * phy.space, phy.height  # len(phy.tips) * 14 + 80
            iface.canvas = toyplot.Canvas(width=w, height=h)
            iface.axes = iface.canvas.cartesian(bounds=(20, w, 0, h - 50 * phy.axis),
                                                ymin=-.5, ymax=iface.tree.ntips - .5)
        elif phy.unro:
            raise NotImplementedError
        else:
            raise RuntimeWarning

        phy.tx += '_TBE' if phy.tbe else '_FBP'

        iface.text = 'plotting %s tree' % phy.tx
        with open(self.wd / 'tree.pickle', 'wb') as out_file:
            pickle.dump(iface.tree, out_file)
        try:
            iface.tup = iface.tree.draw(
                axes=iface.axes, scalebar=phy.axis,
                tip_labels=tl, tip_labels_align=phy.align,
                tip_labels_style={'font-size': '15px'},

                node_style={'stroke': None},
                node_sizes=ns if phy.supp else None,
                node_colors=nc,
                node_labels='support' if phy.supp else None,
                node_labels_style={'font-size': '7px', 'fill': '#FFFFFF',
                                   'baseline-shift': '-1px',
                                   'font-weight': 'bold'},

                tip_labels_colors=tlc,
                layout=ly,

                node_hover=phy.html)

            iface.canvas.style['background-color'] = 'white'

            if phy.circ:
                iface.tup[1].show = False
            elif phy.rect:
                iface.tup[1].y.show = False
                iface.tup[1].x.show = phy.axis
                # iface.tup[1].x.spine.position = 'high' # and swap out the 0,h-50 bounds
                th = iface.tree.treenode.height
                iface.tup[1].x.domain.max = (w - 1200) / w  # max(th / 6, phy.spec * th / 2, phy.align * th / 1.6)
                # 1/5 is not enough for aligned labels
                iface.tup[1].x.domain.min = -th
                # 0 is the right-most tip (not tip label) of the tree
                # therefore this inverted x domain.
            iface.i += 1

            try:
                self.save_tree()
            except ValueError:
                LOG.error('canvas to narrow, widening ...')
                iface.spacespin.props.adjustment.props.value += .1

        except ET.ParseError as ex:
            e = 'XML ParseError. Invalid character in a sample ID? Please check metadata.tsv'
            LOG.error(e)
            GObject.idle_add(self.reset_tree)
            sleep(.1)
            return True
        except (ValueError, IndexError) as ex:
            LOG.exception(ex)
            GObject.idle_add(self.reset_tree)
            sleep(.1)
            return True

        GObject.idle_add(self._do_phy4, errors)
        sleep(.1)
        return True

    def _do_phy4(self, errors):
        """Bounce the plotting to the main loop"""
        data = self.data
        iface = self.iface
        phy = data.phy
        iface.thread.join()
        iface.text = 'placing tree'
        LOG.debug(iface.text)
        self.update(PAGE)
        self.win.show_all()
        if errors:
            self.show_notification('Errors during tree plotting', errors)
            return
        sleep(.1)

        self.load_image(iface.zoomer, PAGE, iface.tree_eventbox,
                        self.wd / (phy.tx + '.png'), h=phy.height)
        self.win.show_all()
        sleep(.1)
        iface.i += 1

        # re-direct to thread
        iface.text = 'plotting MSA'
        iface.thread = threading.Thread(target=self._do_phy5)
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()

    def _do_phy5(self, popgen_markup=None):
        """Plot the MSA matrix"""
        data = self.data
        iface = self.iface
        phy = data.phy

        # the tips are saved top-down, but in the tree they're bottom-up. so reverse
        full_len = sum(phy.g_lens.values()) \
                   + (len(data.genes) - 1) * len(repo.SEP)
        for i, seq in enumerate(phy.seqdata):
            if not seq:
                phy.seqdata[i] = [repo.toint('?')] * full_len

        phy.array = np.array(phy.seqdata[::-1])
        # make gaps transparent
        array = np.ma.masked_where(phy.array > repo.toint('else'), phy.array)

        if popgen_markup:
            # make selection non-transparent
            array = np.flipud(array)
            rows, drop_cols, repeats = popgen_markup
            tree_mask = np.full(array.shape, 1.0)
            # mask the unselected rows
            tree_mask[[i for i in range(array.shape[0])
                       if i not in rows], :] = repo.ALPHA
            # mask the invalid columns
            tree_mask[:, sorted(drop_cols)] = repo.ALPHA

        # adjust maximum size
        scale = 6
        while max(array.shape) * scale > 2 ** 14:
            scale -= 1
        LOG.debug('scaling qal with %d' % scale)

        f = Figure()
        f.set_facecolor('none')
        f.set_figheight(array.shape[0] / repo.DPI * 5)
        f.set_figwidth(array.shape[1] / repo.DPI)

        # this line is now redundant
        f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

        # leave 50px at the bottom for the tree axis and gene legend bar
        b = (50 * phy.axis) / phy.height

        ax = f.add_axes([0, b, 1, 1 - b])
        if popgen_markup:
            mat = ax.matshow(array, alpha=tree_mask, aspect='auto',
                             cmap=ListedColormap(repo.colors),
                             vmin=-.5, vmax=len(repo.colors) - .5)
        else:
            mat = ax.pcolormesh(array, alpha=1, cmap=ListedColormap(repo.colors),
                                vmin=-.5, vmax=len(repo.colors) - .5)
        ax.axis('off')

        # plot a gene marker bar
        if phy.axis and len(data.genes) > 1:
            LOG.debug('adding gene marker bar')
            # build matrix of gene indicators
            gm = [[data.genes.index(g)] * phy.g_lens[g] for g in data.genes]
            # add the spacer
            for i in range(len(data.genes) - 1):
                gm[i] += [-1] * len(repo.SEP)
            gm = [gi for gl in gm for gi in gl]
            # turn into np array
            gm = np.vstack([np.array(gm)] * 2)
            # make spacer transparent
            gm = np.ma.masked_where(gm < 0, gm)
            gene_colors = get_cmap('GnBu', len(data.genes))

            # plot the marker bar onto the MSA graphic
            bax = f.add_axes([0, b * .4, 1, b / 3])
            bar = bax.pcolormesh(gm, cmap=gene_colors)
            # bax.axis('off')
            [bax.spines[t].set_visible(False) for t in ['left', 'right', 'top', 'bottom']]
            bax.yaxis.set_visible(False)
            bax.xaxis.set_ticks_position('bottom')
            bax.tick_params(colors=iface.FG, pad=.2, length=0, labelsize=1)
            bax.xaxis.set_ticks([i for i in range(1, gm.shape[1] - 1) if not i % 100])

            # plot a legend for the gene marker bar
            with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
                iface.text = 'gene marker bar'
                LOG.debug(iface.text)
                Path.mkdir(self.wd / repo.PATHS.phylo_msa.parent, exist_ok=True)
                fig = plt.figure(figsize=(4, .2))
                cax = fig.add_subplot(111)
                cbar = ColorbarBase(
                    ax=cax, cmap=gene_colors, orientation='horizontal',
                    ticks=[(.5 / len(data.genes) + j * 1 / len(data.genes))
                           for j in range(len(data.genes))])
                cbar.ax.set_xticklabels(data.genes)
                fig.savefig(self.wd / repo.PATHS.gbar, transparent=True,
                            bbox_inches='tight', pad_inches=0, dpi=600)
                plt.close(fig)
                del fig, cbar

        iface.i += 1
        iface.fig = f
        iface.scale = scale
        self.save_msa(popgen_markup)

        phy.shape = array.shape[1], phy.height

        if iface.rasterize.props.active:
            iface.text = 'place PNG MSA: %d:%d' % array.shape
            LOG.debug(iface.text)
            sleep(.05)
            self.load_image(iface.zoomer, PAGE, iface.msa_eventbox,
                            self.wd / repo.PATHS.phylo_msa,
                            w=phy.shape[0] * self.get_hadj() * 2, h=phy.shape[1])
        else:
            iface.text = 'place vector'
            LOG.debug(iface.text)
            canvas = FigureCanvas(f)  # a Gtk.DrawingArea
            canvas.set_size_request(max(phy.shape[0] * self.get_hadj() * 2, -1),
                                    phy.shape[1])  # height
            try:
                ch = iface.msa_eventbox.get_child()
                if ch:
                    iface.msa_eventbox.remove(ch)
                iface.msa_eventbox.add(canvas)
            except Exception as ex:
                LOG.error(ex)
        plt.close(f)
        del f
        iface.i += 1
        self.load_colorbar(iface.palplot1)

        if phy.axis and len(data.genes) > 1:
            self.load_colorbar(iface.gbar, gbar=True)
            iface.gbar.set_visible(True)
        else:
            iface.gbar.set_visible(False)

        iface.text = 'idle'
        sleep(.1)
        GObject.idle_add(self._stop_phy)
        return True

    def _stop_phy(self, error=None):
        self.iface.thread.join()
        if error:
            LOG.error(error)
            self.show_notification(error, secs=10)
            self.win.show_all()
            return False
        LOG.info('plotting thread idle')
        self.iface.view_msa_ids.set_model(self.data.tree_anno_model)
        self.win.show_all()
        self.save(silent=True)
        self.data.change_indicator[PAGE] = False
        # *Now* go back to where you came from
        if self.iface.run_after:
            func, args, *_ = *self.iface.run_after, (), ()
            func(*args)
            self.iface.run_after = list()
        return False

    def _on_right_click(self, widget, event, menu):
        # Check if right mouse button was pressed
        if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
            menu.popup(None, None, None, None, event.button, event.time)
            return True  # event has been handled

    def _size_scrollbars(self, pane, pos, scl, scr, sp, tv):
        """
        Set a size request for the scrollbar below the narrower
        side of the GtkPaned and set the other one to hexpand.
        """
        p = pane.get_position()
        w = pane.get_allocated_width()
        p = min(p, w)
        self.iface.ratio = p / w
        if p > w / 2:
            scl.set_size_request(-1, -1)
            scl.set_hexpand(True)
            scr.set_size_request(max(w - p - tv.get_allocated_width(), -1), -1)
            scr.set_hexpand(False)
        elif p < w / 2:
            scl.set_size_request(p, -1)
            scl.set_hexpand(False)
            scr.set_size_request(-1, -1)
            scr.set_hexpand(True)
        sp.set_size_request(tv.get_allocated_width(), -1)

    def _size_pane(self, win, state):
        iface = self.iface
        with iface.tree_pane.handler_block(iface.pane_handler):
            try:
                sleep(.1)
                if win.props.is_maximized:
                    s = win.get_screen()
                    m = s.get_monitor_at_window(s.get_active_window())
                    monitor = s.get_monitor_geometry(m)
                    self.iface.tree_pane.set_position(self.iface.ratio * monitor.width)
                else:
                    self.iface.tree_pane.set_position(self.iface.ratio * win.get_size()[0])
            except Exception as ex:
                LOG.error(ex)

    def reload_ui_state(self, *args):
        data = self.data
        iface = self.iface
        phy = data.phy

        for w_name in ['gap_share', 'unk_share']:
            iface.__getattribute__(w_name).get_adjustment().props.value = \
                phy.__getattribute__(w_name)

        for w_name in ['query', 'exclude']:
            iface.__getattribute__(w_name).set_text(
                phy.__getattribute__(w_name))

        # load plot settings
        [iface.__getattribute__(w_name).set_active(
            phy.__getattribute__(w_name)) for w_name in
            ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
             'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]

        iface.flipspin.props.adjustment.props.value = phy.flip
        iface.distspin.props.adjustment.props.value = phy.dist
        iface.spacespin.props.adjustment.props.value = phy.space
        if data.genes:
            self.init_sp_genes(iface.sp_genes, data.genes, phy.sel_gene)


def _h(a):
    # harmonic number
    b = 0
    for c in range(1, a):
        b += 1 / c
    return b


def _qh(a):
    # quadratic harmonic number
    b = 0
    for c in range(1, a):
        b += 1 / np.power(c, 2)
    return b


def _per_gene_diversity(gene, phy, array, _range):
    conserved, biallelic, polyallelic, unknown, gaps = 0, 0, 0, 0, 0
    start, end = min(_range), max(_range)
    nrows = array.shape[0]
    too_unknown = phy.unk_share * nrows
    too_gappy = phy.gap_share * nrows
    poly = True
    drop = set()

    # iterate over MSA columns, search for segregating sites
    for j in _range:
        if j >= array.shape[1]:
            break
        col = array[:, j]
        if len(set(col)) == 1:
            # only one char at site
            it = col[0]
            if it > 4:
                # all-unknown-site. can happen as we selected a subset of sequences
                unknown += 1
                drop.add(j)
            elif it == 4:
                # all-gap-site
                gaps += 1
                drop.add(j)
            else:
                conserved += 1
        else:  # more than one character in the whole column
            if max(col) > 4 and len(col[col > 4]) > too_unknown:
                # there are too many unknown chars in the column
                unknown += 1
                drop.add(j)
            elif max(col) == 4 and len(col[col >= 4]) > too_gappy:
                # there are too many gaps at the site
                gaps += 1
                drop.add(j)
            else:
                # if gaps or unknown characters are at the site,
                # replace them with the most common nucleotide
                try:
                    col[col >= 4] = np.bincount(col[col < 4]).argmax()
                except ValueError:
                    # there are no nucleotides
                    drop.add(j)
                    continue

                # ordinary site treatment
                num_diff_chars = len(set(col))
                if num_diff_chars > 2:
                    polyallelic += 1
                    if not poly:
                        drop.add(j)
                elif num_diff_chars == 2:
                    biallelic += 1
                else:
                    conserved += 1

    seg_sites = biallelic
    n_sites = end - start + 1 - unknown - gaps
    if poly:
        seg_sites += polyallelic
    else:
        n_sites -= polyallelic

    # crop to allowed sites
    crop = array[:, [i for i in _range if i not in drop]]

    k = 0
    # count pairwise differences
    for combi in itertools.combinations(range(nrows), 2):
        k += np.sum(crop[combi[0]] != crop[combi[1]])

    # binomial coefficient, k = k^ in formula (10) from Tajima1989
    k = k * 2 / nrows / (nrows - 1)
    pi = k / n_sites

    # look up or compute
    if nrows in phy.h:
        a1 = phy.h[nrows]
        a2 = phy.qh[nrows]
    else:
        a1 = _h(nrows)
        phy.h[nrows] = a1
        a2 = _qh(nrows)
        phy.qh[nrows] = a2

    if n_sites == 0:
        return ([gene, n_sites, seg_sites, k, pi, pi, pi,
                 -1, gaps, unknown], drop)

    theta_w = seg_sites / a1 / n_sites  # per site from Yang2014
    unq, idx, count = np.unique(crop, axis=0, return_counts=True, return_index=True)
    repeats = {ri: np.argwhere(np.all(crop == rg, axis=1)).ravel().tolist()
               for rg, ri in zip(unq[count > 1], idx[count > 1])}

    # Tajima's D Formula:
    # d = k - seg_sites/a1 = k - theta_W
    # D = d / sqrt(Var(d))
    try:
        e1 = 1 / a1 * ((nrows + 1) / (3 * nrows - 3) - 1 / a1)
        b2 = (2 * (nrows * nrows + nrows + 3)) / (9 * nrows * (nrows - 1))
        c2 = b2 - ((nrows + 2) / (nrows * a1)) + (a2 / (a1 * a1))
        e2 = c2 / (a1 * a1 + a2)
        tajima = (k - theta_w) / np.sqrt(e1 * seg_sites + e2 * seg_sites * (seg_sites - 1))
    except (ZeroDivisionError, RuntimeWarning) as z:
        tajima = float('inf')

    return ([gene, n_sites, seg_sites, k, pi, theta_w, tajima,
             len(unq), gaps, unknown], drop, repeats)
