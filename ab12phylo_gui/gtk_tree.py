# 2020 Leo Kaindl

import json
import logging
import threading
import toytree
import typing
import toyplot
import numpy as np
from numpy import nan
from pathlib import Path
import xml.etree.ElementTree as ET
from time import sleep
import pandas as pd
from math import ceil
from copy import deepcopy
from toyplot import html, pdf, png, svg, color, data, locator

import gi
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo_gui import static, shared
from ab12phylo_gui.static import PATHS as p, SEP, tohex, \
    blue, red, dark_red, rocket, black

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 7


def init(gui):
    """Initialize the page. Connect buttons"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    # handle edits to the popgen input immediately
    for w_name in ['gap_share', 'unk_share']:
        iface.__getattribute__(w_name).get_adjustment().connect(
            'value-changed', lambda adj: phy.__setattr__(w_name, adj.props.value))

    for w_name in ['query', 'exclude']:
        iface.__getattribute__(w_name).connect(
            'changed', lambda adj: phy.__setattr__(w_name, adj.props.value))

    for w_name in ['matches', 'subtree']:
        iface.__getattribute__(w_name).connect('clicked', start_popgen, gui)

    iface.save_plot.connect('clicked', on_save_plot, gui)

    iface.view_msa_ids.set_model(data.tree_anno_model)
    col = Gtk.TreeViewColumn(title='id', cell_renderer=Gtk.CellRendererText(), text=0)
    col.set_resizable(True)
    iface.view_msa_ids.append_column(col)
    iface.view_msa_ids.set_tooltip_column(1)

    iface.popgen.set_model(data.pop_model)
    set_popgen_columns(iface.popgen)

    # set some plot_menu MenuButton images
    iface.flip.set_image(iface.flipim)
    iface.dist.set_image(iface.distim)
    iface.anno.set_image(iface.annoim)

    # connect zooming and selecting
    sel = iface.view_msa_ids.get_selection()
    sel.set_mode(Gtk.SelectionMode.MULTIPLE)
    sel.connect('changed', shared.keep_visible,
                iface.parallel_tree.props.vadjustment.props, iface.tempspace)
    iface.msa_eventbox.connect('button_press_event', shared.select_seqs, PAGE, iface.zoom,
                               iface.view_msa_ids, iface.tempspace)  # in-preview selection
    iface.msa_eventbox.connect('scroll-event', shared.xy_scale, gui, PAGE)  # zooming

    iface.tree = None
    iface.phy_seen = False


def start_popgen(widget, gui):
    """Prepare the plotting Thread"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return
    phy.mode = widget.get_name()
    LOG.debug('start_popgen', phy.mode)

    data.pop_model.clear()

    # TODO
    set_popgen_columns(iface.popgen)
    return


def set_popgen_columns(tv):
    if len(tv.get_model()) == 0 or tv.get_n_columns > 0:
        return
    for i, ti in enumerate(['gene', 'sites', 'S', 'k', 'π', 'θ',
                            'Tajima\'s D', 'unique', 'gap', 'unknown']):
        tv.append_column(Gtk.TreeViewColumn(
            title=ti, cell_renderer=Gtk.CellRendererText(), text=i))


def on_save_plot(wi, gui):
    """Save the current plot in different formats than before"""
    LOG.debug('on_save_plot')


def start_phy(gui, run_after=None):
    """Get settings and block GUI"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return
    LOG.debug('start_plotting')

    # parse plot settings only when needed
    [phy.__setattr__(w_name, iface.__getattribute__(w_name).get_active()) for w_name in
     ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec', 'axis', 'pmsa', 'pdf', 'svg', 'png', 'nwk']]
    phy.flip = iface.flipspin.props.adjustment.props.value
    phy.dist = iface.distspin.props.adjustment.props.value
    phy.sel_gene = iface.sp_genes.get_active_text()

    # re-direct to thread
    iface.thread = threading.Thread(target=do_phy1, args=[gui])
    iface.run_after = run_after
    iface.k = (phy.rect + phy.circ + phy.unro) * 2 + 6
    iface.i = 0
    iface.text = 'read tree'
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    return
    # return to main loop


def do_phy1(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    tree_file = gui.wd / p.fbp if phy.fbp else gui.wd / p.tbe
    phy.scale = .1 if phy.fbp else 10

    # read tree file
    iface.tree = toytree.tree(open(tree_file, 'r').read(), tree_format=0)

    # fetch tips labels
    phy.tips = list(reversed(iface.tree.get_tip_labels()))
    iface.i += 1

    iface.text = 'read metadata'
    # read in tabular data
    if 'df' not in iface.tempspace:
        iface.tempspace.df = pd.read_csv(gui.wd / p.tsv, sep='\t', dtype={'id': str})
        iface.tempspace.df.set_index('id', inplace=True)
    df = iface.tempspace.df.copy(deep=True)

    r, b, e = 'reference_species', 'BLAST_species', 'extra_species'
    if r in df:
        # overwrite BLAST with reference
        df.loc[df[r].notna(), b] = df.loc[df[r].notna(), r]
    if b in df:
        df = df.rename(columns={b: 'species'})
        phy.did_BLAST = True
    # crop to relevant columns
    df = df[[c for c in ['gene', 'accession', 'species', 'pid'] if c in df]]

    df = df.fillna(value={'pid': 0})
    pid = 'pid' in df
    if phy.sel_gene != 'best phred' or len(data.genes) == 1:
        df = df.loc[df['gene'] == phy.sel_gene]
        df = df.loc[phy.tips]
        phy.ndf = df.to_dict(orient='index')
    else:
        # split into per-gene dataframes
        dfs = [i[1] for i in df.groupby(df['gene'])]
        # find row with highest pid from all genes
        phy.ndf = dict()
        for t in phy.tips:  # TODO handle case without pid and what about references?
            rows = [d.loc[t] for d in dfs]
            if pid:
                phy.ndf[t] = [dict(r) for r in rows if r.pid == max([q.pid for q in rows])][0]
            else:
                phy.ndf[t] = [dict(r) for r in rows][0]

    iface.text = 'fill id column'
    data.tree_anno_model.clear()
    iface.i += 1
    sleep(.1)
    GObject.idle_add(do_phy2, gui)
    return True


def do_phy2(gui):
    """Bounce to idle to size the id column"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    # write samples to model
    for t in phy.tips:
        sp = phy.ndf[t].get('species', '')
        if not sp:
            pass
        elif pd.isna(sp):
            sp = ''
            phy.ndf[t].pop('species')
            phy.ndf[t].pop('pid')
        # pid = ndf[t]['pid']
        # if pd.isna(pid) or pid == 0:
        #     ndf[t].pop('pid')
        c = iface.FG if 'accession' not in phy.ndf[t] else iface.BLUE
        data.tree_anno_model.append([t, sp, c, iface.BG])

    gui.win.show_all()
    sleep(.1)

    # re-direct to thread
    iface.thread = threading.Thread(target=do_phy3, args=[gui])
    iface.thread.start()


def do_phy3(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.text = 'write new MSA'
    errors = list()

    # read MSA
    records = {r.id: r.seq for r in SeqIO.parse(gui.wd / p.msa, 'fasta')}
    phy.seqdata = list()
    # write re-named re-ordered MSA
    with open(gui.wd / p.msa_anno, 'w') as new_msa:
        for tip in phy.tips:
            # get seq and cut out artificial separator
            seq = str(records[tip]).replace(SEP, '')
            # get metadata
            entry = phy.ndf[tip]

            des = entry.get('species', '')
            pid = '%.2f' % entry['pid'] if 'pid' in entry else ''
            if des and pid:
                des += '; pid %s' % pid

            phy.seqdata.append(records[tip])

            # rename REFs
            if tip.startswith('REF_'):
                tip = entry['accession']

            SeqIO.write(SeqRecord(Seq(seq), id=tip, description=des), new_msa, 'fasta')
    LOG.debug('wrote updated MSA')
    iface.i += 1

    # get lengths of genes
    phy.g_lens = [(gene, len(seq)) for gene, seq in
                  zip(data.genes, records[next(iter(records))].split(SEP))]
    multi = True if len(data.genes) > 1 else False

    species, cat, pid, strain = 'species', 'cat', 'pid', 'strain'

    iface.text = 'iterate tree: 1'
    # rename reference nodes and add features
    for node in iface.tree.treenode.traverse():
        if node.is_leaf():
            entry = phy.ndf[node.name]

            if node.name.startswith('REF_'):
                if multi and strain in entry[species]:
                    node.name = entry[species][entry[species].find(strain) + 7:]
                else:
                    node.name = entry['accession']
                node.add_feature(cat, blue)
                node.add_feature(species, entry[species])
                node.add_feature(pid, 0)  # this is the index of blue in
                # the rocket palette, which is used to display BLAST pid

            elif species in entry:
                node.add_feature(cat, black)  # black
                node.add_feature(species, entry[species])
                node.add_feature(pid, 1 + ceil(min(100 - entry[pid], 20) / 2))  # map to the palette

            else:
                node.add_feature(cat, rocket[1] if not phy.did_BLAST else dark_red)
                node.add_feature(pid, 1 if not phy.did_BLAST else 11)
                # dark gray if no BLAST, else bad light color signifying no hit

            node.add_feature('size', 0)
            node.add_feature('color', blue)  # TODO just so the field exists
        else:
            # use support values for node size
            min_size = 2
            node.add_feature('size', min_size + float(node.support) * phy.scale if node.support else 0)
            node.add_feature('color', blue if node.size > phy.flip + min_size else red)
            node.add_feature(pid, -1)  # just so the field exists
            node.add_feature(cat, 0)  # just so the field exists
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
        tf = gui.wd / p.fbpn if phy.fbp else gui.wd / p.tben
        exp_tree.write(tf, tree_format=0)

    # tree plotting tweaks
    iface.text = 'iterate tree: 2'
    for node in iface.tree.treenode.traverse():
        node.dist = max(node.dist, phy.dist)
        if phy.supp:
            node.support = int(round(node.support * 100))
    iface.i += 1
    # get model height
    phy.height = iface.view_msa_ids.get_allocated_height()
    # plot tree and msa
    if phy.circ:
        iface.text = 'plotting circular'
        try:
            iface.ctup = iface.tree.draw(
                width=800, height=800, scalebar=phy.axis,
                node_sizes=list(iface.tree.get_node_values('size', 1, 1)),
                node_colors=list(iface.tree.get_node_values('color', 1, 1)),
                tip_labels=True, tip_labels_align=False,
                tip_labels_colors=[n for n in list(iface.tree.get_node_values('pid', 1, 1))
                                   if n != -1][::-1],
                layout='c')
            iface.ctup[0].style['background-color'] = 'white'
            iface.ctup[1].show = False
            iface.i += 1

            if phy.png:
                png.render(iface.ctup[0], str(gui.wd / 'circular.png'), scale=1.6)
            if phy.pdf:
                pdf.render(iface.ctup[0], str(gui.wd / 'circular.pdf'))
            if phy.svg:
                svg.render(iface.ctup[0], str(gui.wd / 'circular.svg'))
            iface.i += 1

        except ET.ParseError as ex:
            e = 'XML ParseError. Invalid character in a sample ID? Please check metadata.tsv'
            LOG.error(e)
            errors.append(e)
        except ValueError as ex:
            LOG.exception(ex)
            errors.append(ex)
    elif phy.rect:
        try:
            iface.text = 'plotting rectangular'
            # get dim for canvas
            w, h = 1200, len(phy.tips) * 14 + 80
            if phy.supp:
                a = [color.rgb(*a) for a in
                     list(iface.tree.get_node_values('color', 1, 1))]

                iface.rtup = iface.tree.draw(width=w, height=h, scalebar=phy.axis, tip_labels=True,
                                             node_labels='support',
                                             node_labels_style={'font-size': '6px', 'fill': '#FFFFFF',
                                                                'baseline-shift': '-1px',
                                                                'font-weight': 'bold'},
                                             node_sizes=list(iface.tree.get_node_values('size', 1, 1)),
                                             node_colors=[color.rgb(*a[:-1]) for a in
                                                          list(iface.tree.get_node_values('color', 1, 1))])
                # tip_labels_colors=[rocket[n] for n in list(
                #     iface.tree.get_node_values('pid', 1, 1)) if n != -1][::-1])
            else:
                print(toytree.__version__)
                a = [color.css(n) for n in list(
                    iface.tree.get_node_values('pid', 1, 1)) if n != -1][::-1]

                iface.rtup = iface.tree.draw(width=w, height=h, scalebar=phy.axis, tip_labels=True,
                                             node_sizes=list(iface.tree.get_node_values('size', 1, 1)),
                                             node_style={'stroke': None},
                                             node_colors=[color.css(a) for a in
                                                          list(iface.tree.get_node_values('color', 1, 1))])  #
                # tip_labels_colors=[color.rgb(rocket[n][0], rocket[n][1], rocket[n][2])
                #                    for n in list(
                #         iface.tree.get_node_values('pid', 1, 1)) if n != -1][::-1])
            iface.rtup[0].style['background-color'] = 'white'
            iface.rtup[1].y.show = False
            iface.rtup[1].x.show = True
            iface.rtup[1].x.domain.max = iface.tree.treenode.height / 5  # 0 is right-most tip of tree. -> divide space!
            iface.i += 1

            if phy.png:
                png.render(iface.rtup[0], str(gui.wd / 'rectangular.png'), scale=1.6)
            if phy.pdf:
                pdf.render(iface.rtup[0], str(gui.wd / 'rectangular.pdf'))
            if phy.svg:
                svg.render(iface.rtup[0], str(gui.wd / 'rectangular.svg'))
            iface.i += 1

        except ET.ParseError as ex:
            e = 'XML ParseError. Invalid character in a sample ID? Please check metadata.tsv'
            LOG.error(e)
            errors.append(e)
        except ValueError as ex:
            LOG.exception(ex)
            errors.append(ex)
    elif phy.unro:
        pass
    else:
        raise NotImplementedError

    iface.text = 'idle'
    iface.frac = 1
    sleep(.1)
    GObject.idle_add(refresh, gui)
    return True


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    # algo = static.toalgo(iface.ml_stack.get_visible_child_name())
    # check_MSA(None, gui, algo)
    if not iface.phy_seen:
        phy.query = ''  # TODO delete later
        phy.exclude = ''  # TODO delete later

        reload_ui_state(gui)
        shared.load_colorbar(iface.palplot1, gui.wd)
        init_sp_genes(iface.sp_genes, data.genes, phy.sel_gene)
        start_phy(gui)
        iface.phy_seen = True
        iface.did_BLAST = False  # TODO delete later


def init_sp_genes(sp_genes, genes, sel_gene=None):
    """
    Init the selector for which gene to draw species annotations from
    :param sp_genes: the GtkComboBoxText to display the genes
    :param genes: the iterable of genes in the analysis
    :param sel_gene: an optional selected gene
    :return: sel_gene: if not set before, is 'best phred' now
    """
    sp_genes.remove_all()
    genes = list(genes)
    if len(genes) > 1:
        genes.insert(0, 'best phred')
    [sp_genes.append_text(gene) for gene in genes]
    if sel_gene and sel_gene in genes:
        idx = genes.index(sel_gene)
    else:
        idx = 0
        sel_gene = 'best phred'
    sp_genes.set_active(idx)
    return sel_gene


def reload_ui_state(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    for w_name in ['gap_share', 'unk_share']:
        iface.__getattribute__(w_name).get_adjustment().props.value = \
            phy.__getattribute__(w_name)

    for w_name in ['query', 'exclude']:
        iface.__getattribute__(w_name).set_text(phy.__getattribute__(w_name))

    if data.genes:
        init_sp_genes(iface.sp_genes, data.genes, phy.sel_gene)
    else:
        raise NotImplementedError
