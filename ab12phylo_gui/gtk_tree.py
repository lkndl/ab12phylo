# 2020 Leo Kaindl

import itertools
import logging
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
import re
import svgutils.transform as sg
import toyplot
import toytree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.cm import get_cmap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from toyplot import html, pdf, png, svg, color

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, Gdk

from ab12phylo_gui import shared
from ab12phylo_gui.static import PATHS as p, \
    seqtoint, tohex, toint, SEP, DPI, \
    colors, blue, red, dark_red, rocket, black

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 7
UI_INFO = """
<ui>
  <popup name='context_menu'>
    <menuitem action='root' />
    <menuitem action='drop' />
    <menuitem action='collapse' />
    <menuitem action='extract' />
  </popup>
</ui>
"""


def init(gui):
    """Initialize the page. Connect buttons"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    # context menu
    action_group = Gtk.ActionGroup(name='context_menu')
    action_group.add_actions(
        [('root', None, 'Root', None,
          'Root the tree at the selected node(s).', gui.root_extract),
         ('drop', None, 'Drop', '<Delete>',
          'Drop the selected nodes from the tree.', gui.drop_collapse),
         ('collapse', None, 'Collapse', None,
          'Collapse all descendants to a placeholder.', gui.drop_collapse),
         ('extract', None, 'Extract', None,
          'Drop all other nodes from the tree.', gui.root_extract)])

    uimanager = Gtk.UIManager()
    uimanager.add_ui_from_string(UI_INFO)
    gui.win.add_accel_group(uimanager.get_accel_group())
    uimanager.insert_action_group(action_group)
    iface.contextmenu = uimanager.get_widget('/context_menu')

    for w_name in ['tree_eventbox', 'view_msa_ids', 'msa_eventbox']:
        iface.__getattribute__(w_name).connect(
            'button-press-event', _on_right_click, iface.contextmenu)

    for w_name in ['gap_share', 'unk_share']:
        iface.__getattribute__(w_name).get_adjustment().connect(
            'value-changed', lambda adj: phy.__setattr__(w_name, adj.props.value))

    for w_name in ['query', 'exclude']:
        iface.__getattribute__(w_name).connect('focus_out_event', select_popgen, gui)

    iface.save_plot.connect('clicked', lambda *args: on_save_tree(gui))

    iface.view_msa_ids.set_model(data.tree_anno_model)
    col = Gtk.TreeViewColumn(title='id', text=0, foreground=2,
                             cell_renderer=Gtk.CellRendererText())
    col.set_resizable(True)
    iface.view_msa_ids.append_column(col)
    iface.view_msa_ids.set_tooltip_column(1)

    iface.popgen.set_model(data.pop_model)
    init_popgen_columns(iface.popgen)

    # set some plot_menu MenuButton images
    iface.flip.set_image(iface.flipim)
    iface.dist.set_image(iface.distim)
    iface.anno.set_image(iface.annoim)
    iface.filetypes.set_image(iface.saveim)

    # connect zooming and selecting
    iface.tree_sel = iface.view_msa_ids.get_selection()
    iface.tree_sel.set_mode(Gtk.SelectionMode.MULTIPLE)
    iface.tree_sel.connect('changed', shared.keep_visible,
                           iface.parallel_tree.props.vadjustment.props, iface.tempspace)
    iface.msa_eventbox.connect_after('button_press_event', shared.select_seqs, PAGE, iface.zoomer,
                                     iface.view_msa_ids, iface.tempspace)  # in-preview selection
    iface.tree_eventbox.connect_after('button_press_event', shared.select_seqs, PAGE, iface.zoomer,
                                      iface.view_msa_ids, iface.tempspace)  # in-preview selection
    iface.tree_pane.connect('scroll-event', shared.xy_scale, gui, PAGE)  # zooming
    iface.tree_pane.connect('notify::position', _size_scrollbars,
                            iface.tree_left_scrollbar, iface.tree_right_scrollbar,
                            iface.tree_spacer, iface.view_msa_ids, )

    phy.tx = 'circular' if phy.circ else 'rectangular'
    phy.tx += '_TBE' if phy.tbe else '_FBP'

    iface.tree = None


def _on_right_click(widget, event, menu):
    # Check if right mouse button was pressed
    if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
        menu.popup(None, None, None, None, event.button, event.time)
        return True  # event has been handled


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    reload_ui_state(gui)
    iface.tree_pane.set_position(301)

    if (gui.wd / phy.tx).with_suffix('.png').is_file() \
            and (gui.wd / p.phylo_msa).is_file():
        # load tree PNG
        shared.load_image(iface.zoomer, PAGE, iface.tree_eventbox,
                          gui.wd / (phy.tx + '.png'), h=phy.shape[1])
        # load MSA PNG
        shared.load_image(iface.zoomer, PAGE, iface.msa_eventbox, gui.wd / p.phylo_msa,
                          w=phy.shape[0] * shared.get_hadj(iface), h=phy.shape[1])
        # load colorbars
        shared.load_colorbar(iface.palplot1, gui.wd)
        if phy.axis and len(data.genes) > 1:
            shared.load_colorbar(iface.gbar, gui.wd, gbar=True)
            iface.gbar.set_visible(True)
        gui.win.show_all()
        return

    start_phy(gui)
    gui.win.show_all()


def select_popgen(widget, ev, gui):
    """Select samples for popgen calculations by name / species matching"""
    LOG.debug('select_popgen')
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    phy.__setattr__(widget.get_name(), widget.get_text())

    if not phy.query and not phy.exclude:
        return
    mo = data.tree_anno_model
    rgx = lambda c: re.compile(r'|'.join(
        ['.*' + re.escape(word.strip()) + '.*' for word in
         c.strip().strip(',').split(',')]))
    search_space = [_id + sp for _id, sp in shared.get_column(mo, (0, 1))]
    sel = list()

    # first match
    if not phy.query:
        sel = list(range(len(mo)))
    else:
        sel = [i for i, label in enumerate(search_space) if rgx(phy.query).match(label)]

    # then exclude
    if phy.exclude:
        sel = [i for i, label in enumerate(search_space) if
               i in sel and not rgx(phy.exclude).match(label)]

    # then select in the model
    iface.tree_sel.unselect_all()
    [iface.tree_sel.select_path(i) for i in sel]


def expand_popgen(gui):
    # TODO popgen calc
    # TODO check rooting, dropping, extracting, replacing. and make chainable

    return


def do_popgen(widget, gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    data.pop_model.clear()
    init_popgen_columns(iface.popgen)

    pass


def init_popgen_columns(tv):
    if len(tv.get_model()) == 0 or tv.get_n_columns > 0:
        return
    for i, ti in enumerate(['gene', 'sites', 'S', 'k', 'π', 'θ',
                            'Tajima\'s D', 'unique', 'gap', 'unknown']):
        tv.append_column(Gtk.TreeViewColumn(
            title=ti, cell_renderer=Gtk.CellRendererText(), text=i))


def init_sp_genes(sp_genes, genes, sel_gene=None):
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


def on_save_tree(gui):
    """Save the current tree in different formats"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    LOG.debug('rendering tree')
    if 'canvas' not in iface:
        start_phy(gui)
        return
    tx = gui.wd / phy.tx
    iface.text = 'rendering png'
    png.render(iface.canvas, str(tx.with_suffix('.png')), scale=1.6)
    if phy.pdf:
        iface.text = 'rendering pdf'
        pdf.render(iface.canvas, str(tx.with_suffix('.pdf')))
    if phy.svg or phy.pmsa:
        iface.text = 'rendering svg'
        svg.render(iface.canvas, str(tx.with_suffix('.svg')))
    if phy.html:
        iface.text = 'rendering html'
        html.render(iface.canvas, str(tx.with_suffix('.html')))
    return


def on_save_msa(gui):
    """Save the current msa in different formats"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    LOG.debug('rendering msa')
    Path.mkdir(gui.wd / p.phylo_msa.parent, exist_ok=True)
    msa = gui.wd / p.phylo_msa.parent / p.phylo_msa.stem
    dpi = iface.scale * DPI
    iface.text = 'rendering png'
    iface.fig.savefig(msa.with_suffix('.png'), transparent=True, dpi=dpi)
    if phy.svg:
        iface.text = 'rendering svg'
        iface.fig.savefig(msa.with_suffix('.svg'), transparent=True)
    if phy.pdf or phy.pmsa:
        iface.text = 'rendering pdf'
        iface.fig.savefig(msa.with_suffix('.pdf'), transparent=True, dpi=dpi)

    if phy.pmsa:
        iface.text = 'collating tree and msa'
        tree = sg.fromfile((gui.wd / phy.tx).with_suffix('.svg'))
        m = sg.fromfile(msa.with_suffix('.svg'))
        mw, mh = [float(i[:-2]) * 96 / 72 for i in m.get_size()]  # now in px
        pt = tree.getroot()
        tw, th = [float(i[:-2]) for i in tree.get_size()]  # already in px
        pm = m.getroot()
        pm.moveto(tw, 0, th / mh / 6, th / mh * 96 / 72)
        fig = sg.SVGFigure()  # '%.0fpx' % (tw + th / mh / 8), '%.0fpx' % th)
        fig.append([pt, pm])
        fig.save(gui.wd / (phy.tx + '_msa.svg'))


def start_phy(gui, kwargs=None):
    """Get settings and block GUI"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return
    LOG.debug('start_plotting')
    if not kwargs:
        kwargs = dict()

    # parse plot settings only when needed
    [phy.__setattr__(w_name, iface.__getattribute__(w_name).get_active()) for w_name in
     ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
      'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]
    phy.flip = iface.flipspin.props.adjustment.props.value
    phy.dist = iface.distspin.props.adjustment.props.value
    phy.sel_gene = iface.sp_genes.get_active_text()

    # re-direct to thread
    iface.thread = threading.Thread(target=do_phy1, args=[gui, kwargs])
    iface.k = (phy.rect + phy.circ + phy.unro) * 2 + 9
    iface.i = 0
    iface.text = 'read tree'
    sleep(.1)
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    return
    # return to main loop


def do_phy1(gui, kwargs):
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    tree_file = gui.wd / p.fbp if phy.fbp else gui.wd / p.tbe

    # read tree file
    iface.tree = toytree.tree(open(tree_file, 'r').read(), tree_format=0)

    if 'collapse' in kwargs:
        for idx in kwargs['collapse']:
            iface.tree.idx_dict[idx].add_sister(name='%d_replaced' % idx, dist=1)
            iface.tree.idx_dict[idx].detach()
        # saving the node names now saves some work
        kwargs['collapse'] = ['%d_replaced' % idx for idx in kwargs['collapse']]
    else:
        # save empty list
        kwargs['collapse'] = list()

    if 'drop' in kwargs:
        [iface.tree.idx_dict[idx].detach() for idx in kwargs['drop']]

    # outgroup rooting
    if 'root' in kwargs:
        iface.tree = iface.tree.root(names=kwargs['root'])

    if 'extract' in kwargs:
        treenode = iface.tree.idx_dict[iface.tree.
            get_mrca_idx_from_tip_labels(names=kwargs['extract'])]
        iface.tree = toytree.tree(treenode.write())

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
    if phy.sel_gene != 'best pid' or len(data.genes) == 1:
        df = df.loc[df['gene'] == phy.sel_gene]
        df = df.loc[phy.tips]
        phy.ndf = df.to_dict(orient='index')
    else:
        # split into per-gene dataframes
        dfs = [i[1] for i in df.groupby(df['gene'])]
        # find row with highest pid from all genes
        phy.ndf = dict()
        for t in phy.tips:
            rows = [d.loc[t] for d in dfs]
            if pid:
                phy.ndf[t] = [dict(r) for r in rows if r.pid == max([q.pid for q in rows])][0]
            else:
                phy.ndf[t] = [dict(r) for r in rows][0]

    iface.i += 1
    GObject.idle_add(do_phy2, gui)
    sleep(.1)
    return True


def do_phy2(gui):
    """Bounce to idle to size the id column"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.thread.join()
    iface.text = 'fill id column'
    LOG.debug(iface.text)
    data.tree_anno_model.clear()

    # write samples to model
    for t in phy.tips:
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
            c = iface.BLUE
            sp = t + ': ' + sp
            t = nd['accession']
        else:
            c = None  # iface.FG
        data.tree_anno_model.append([t, sp, c, None])

    for wi in [iface.tree_eventbox, iface.msa_eventbox]:
        wi.get_child().destroy()
    iface.view_msa_ids.realize()
    sleep(.05)

    # re-direct to thread
    iface.thread = threading.Thread(target=do_phy3, args=[gui])
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()


def do_phy3(gui):
    """Write new, annotated MSA, annotated nwk tree and plot tree"""
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
            pid = ('%.2f' % entry['pid']).rstrip('0').rstrip('.')
            if des and pid != '0':
                des += '; pid %s' % pid

            phy.seqdata.append(seqtoint(records[tip]))

            # rename REFs
            if tip.startswith('REF_'):
                tip = entry['accession']

            SeqIO.write(SeqRecord(Seq(seq), id=tip, description=des), new_msa, 'fasta')
    LOG.debug('wrote updated MSA')
    iface.i += 1

    # get lengths of genes
    phy.g_lens = {gene: len(seq) for gene, seq in
                  zip(data.genes, records[next(iter(records))].split(SEP))}
    multi = True if len(data.genes) > 1 else False

    # cat: REF_ -> 1, no_BLAST_hit -> -1, regular -> 0
    # pid: discrete index in rocket palette
    # confident: 1 if >= flip threshold -> blue, else 0 -> red

    iface.text = 'iterate tree: 1'
    cat_ignore = -1
    # rename reference nodes and add features
    for node in iface.tree.treenode.traverse():
        if not node.support:
            node.add_feature('support', 0)
        elif type(node.support) != float:
            node.support = float(node.support)
        if node.is_leaf():
            entry = phy.ndf[node.name]

            if node.name.startswith('REF_'):
                if multi and 'strain' in entry['species']:
                    strain = entry['species'][entry['species'].find('strain') + 7:]
                    phy.ndf['strain'] = phy.ndf.pop(node.name)
                    node.name = strain
                else:
                    acc = entry['accession']
                    phy.ndf[acc] = phy.ndf.pop(node.name)
                    node.name = acc
                node.add_feature('cat', 1)  # blue
                node.add_feature('species', entry['species'])
                node.add_feature('pid', 0)  # this is the index of blue in
                # the rocket palette, which is used to display BLAST pid

            elif 'species' in entry:
                node.add_feature('species', entry['species'])
                node.add_feature('cat', 0 if entry['species'].strip() else 2)  # black else dark_red
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
            if phy.fbp:
                node.support /= 100
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
        tf = gui.wd / p.fbpn if phy.fbp else gui.wd / p.tben
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
    # get model height
    phy.height = iface.view_msa_ids.get_allocated_height() + 50 * phy.axis

    # pre-calculate plot attributes
    # internal node size -> support
    min_size = .2  # TODO fails on drop
    ns = [min_size + s / 10 for s in iface.tree.get_node_values('support', 1, 1)]

    # internal node color -> support
    red_blue = [color.rgb(*c) for c in [red, blue]]
    nc = [red_blue[i] for i in iface.tree.get_node_values('confident', 1, 1)]

    # primary tip label color -> category
    cat_rgb = [color.rgb(*c) for c in [black, blue, dark_red]]
    tlc = [cat_rgb[n] for n in list(iface.tree.get_node_values(
        'cat', 1, 1)) if n != cat_ignore][::-1]

    # secondary tip label color -> BLAST pid on species
    tl = iface.tree.get_tip_labels()
    if phy.spec:
        # rocket_rgb = [color.rgb(*c) for c in rocket]
        rocket_css = [tohex(c) for c in rocket]
        tl = ['%s <i><span style="fill:%s">%s</span></i>'
              % (t.name, rocket_css[t.pid], t.species) for t in
              iface.tree.treenode.traverse('postorder') if t.is_leaf()][::-1]

    # plot tree and msa
    if phy.circ:
        phy.tx = 'circular'
        ly = 'c'
        w, h = 800, 800
        iface.canvas = toyplot.Canvas(width=w, height=h)
        iface.axes = iface.canvas.cartesian(bounds=(0, w, 0, h))
    elif phy.rect:
        phy.tx = 'rectangular'
        ly = 'r'
        # get dim for canvas
        w, h = 1200, phy.height  # len(phy.tips) * 14 + 80
        iface.canvas = toyplot.Canvas(width=w, height=h)
        iface.axes = iface.canvas.cartesian(bounds=(20, w, 0, h - 50 * phy.axis),
                                            ymin=-.5, ymax=iface.tree.ntips - .5)
    elif phy.unro:
        raise NotImplementedError
    else:
        assert False

    phy.tx += '_TBE' if phy.tbe else '_FBP'

    iface.text = 'plotting %s tree' + phy.tx
    try:
        iface.tup = iface.tree.draw(
            axes=iface.axes, scalebar=phy.axis,
            tip_labels=tl, tip_labels_align=phy.align,

            node_style={'stroke': None},
            node_sizes=ns if phy.supp else None,
            node_colors=nc,
            node_labels='support' if phy.supp else None,
            node_labels_style={'font-size': '6px', 'fill': '#FFFFFF',
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
            iface.tup[1].x.domain.max = th / 5 + phy.align * th / 3
            # 1/5 is not enough for aligned labels
            iface.tup[1].x.domain.min = -th
            # 0 is the right-most tip (not tip label) of the tree
            # therefore this inverted x domain.
        iface.i += 1

        on_save_tree(gui)
        iface.i += 1

    except ET.ParseError as ex:
        e = 'XML ParseError. Invalid character in a sample ID? Please check metadata.tsv'
        LOG.error(e)
        errors.append(e)
    except ValueError as ex:
        LOG.exception(ex)
        errors.append(ex)

    GObject.idle_add(do_phy4, gui, errors)
    sleep(.1)
    return True


def do_phy4(gui, errors):
    """Bounce the plotting to the main loop"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.thread.join()
    iface.text = 'placing tree'
    LOG.debug(iface.text)
    shared.update(iface, PAGE)
    gui.win.show_all()
    if errors:
        shared.show_notification(gui, 'Errors during tree plotting', errors)
        return
    sleep(.1)

    shared.load_image(iface.zoomer, PAGE, iface.tree_eventbox,
                      gui.wd / (phy.tx + '.png'), h=phy.height)
    gui.win.show_all()
    sleep(.1)
    iface.i += 1

    # re-direct to thread
    iface.text = 'plotting MSA'
    iface.thread = threading.Thread(target=do_phy5, args=[gui])
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()


def do_phy5(gui):
    """Plot the MSA matrix"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    array = np.array(phy.seqdata)
    # make gaps transparent
    array = np.ma.masked_where(array > toint('else'), array)

    # adjust maximum size
    scale = 6
    while max(array.shape) * scale > 2 ** 14:
        scale -= 1
    LOG.debug('scaling qal with %d' % scale)

    f = Figure()
    f.set_facecolor('none')
    f.set_figheight(array.shape[0] / DPI)
    f.set_figwidth(array.shape[1] / DPI)

    # this line is now redundant
    f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    # leave 50px at the bottom for the tree axis and gene legend bar
    b = (50 * phy.axis) / phy.height

    ax = f.add_axes([0, b, 1, 1 - b])
    mat = ax.pcolormesh(array, alpha=1, cmap=ListedColormap(colors),
                        vmin=-.5, vmax=len(colors) - .5)  # , aspect='auto')
    ax.axis('off')

    # plot a gene marker bar
    if phy.axis and len(data.genes) > 1:
        LOG.debug('adding gene marker bar')
        # print(phy.g_lens)
        # build matrix of gene indicators
        gm = [[data.genes.index(g)] * phy.g_lens[g] for g in data.genes]
        # add the spacer
        for i in range(len(data.genes) - 1):
            gm[i] += [-1] * len(SEP)
        gm = [gi for gl in gm for gi in gl]
        # turn into np array
        gm = np.vstack([np.array(gm)] * 2)
        # make spacer transparent
        gm = np.ma.masked_where(gm < 0, gm)
        gene_colors = get_cmap('GnBu', len(data.genes))

        # plot the marker bar onto the MSA graphic
        bax = f.add_axes([0, b / 3, 1, b / 3])
        bar = bax.pcolormesh(gm, cmap=gene_colors)  # , aspect='auto')
        bax.axis('off')

        # plot a legend for the gene marker bar
        with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
            iface.text = 'gene marker bar'
            LOG.debug(iface.text)
            fig = plt.figure(figsize=(4, .2))
            cax = fig.add_subplot(111)
            cbar = ColorbarBase(ax=cax, cmap=gene_colors, orientation='horizontal',
                                ticks=[(.5 / len(data.genes) + j * 1 / len(data.genes))
                                       for j in range(len(data.genes))])
            cbar.ax.set_xticklabels(data.genes)
            fig.savefig(gui.wd / p.gbar, transparent=True,
                        bbox_inches='tight', pad_inches=0, dpi=600)
            plt.close(fig)
            del fig, cbar

    iface.i += 1
    iface.fig = f
    iface.scale = scale
    on_save_msa(gui)

    phy.shape = array.shape[1], phy.height

    if iface.rasterize.props.active:
        iface.text = 'place PNG MSA: %d:%d' % array.shape
        LOG.debug(iface.text)
        sleep(.05)
        shared.load_image(iface.zoomer, PAGE, iface.msa_eventbox, gui.wd / p.phylo_msa,
                          w=phy.shape[0] * shared.get_hadj(iface), h=phy.shape[1])
    else:
        iface.text = 'place vector'
        LOG.debug(iface.text)
        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
        canvas.set_size_request(w=phy.shape[0] * shared.get_hadj(iface), h=phy.shape[1])
        try:
            iface.qal_eventbox.get_child().destroy()
            iface.qal_eventbox.add(canvas)
        except Gtk.Error as ex:
            LOG.error(ex)
    plt.close(f)
    del f
    iface.i += 1
    shared.load_colorbar(iface.palplot1, gui.wd)
    if phy.axis and len(data.genes) > 1:
        shared.load_colorbar(iface.gbar, gui.wd, gbar=True)
        iface.gbar.set_visible(True)

    # save for re-use in popgen part?
    phy.array = array

    iface.text = 'idle'
    sleep(.1)
    GObject.idle_add(stop_phy, gui)
    return True


def stop_phy(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.thread.join()
    LOG.info('phylo thread idle')
    iface.view_msa_ids.set_model(data.tree_anno_model)
    gui.win.show_all()
    data.change_indicator[PAGE] = False
    return False


def _size_scrollbars(pane, pos, scl, scr, sp, tv):
    """
    Set a size request for the scrollbar below the narrower
    side of the GtkPaned and set the other one to hexpand.
    """
    p = pane.get_position()
    w = pane.get_allocated_width()
    if p > w / 2:
        scl.set_size_request(-1, -1)
        scl.set_hexpand(True)
        scr.set_size_request(pane.get_allocated_width() - p - tv.get_allocated_width(), -1)
        scr.set_hexpand(False)
    elif p < w / 2:
        scl.set_size_request(p, -1)
        scl.set_hexpand(False)
        scr.set_size_request(-1, -1)
        scr.set_hexpand(True)
    sp.set_size_request(tv.get_allocated_width(), -1)


def reload_ui_state(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    for w_name in ['gap_share', 'unk_share']:
        iface.__getattribute__(w_name).get_adjustment().props.value = \
            phy.__getattribute__(w_name)

    for w_name in ['query', 'exclude']:
        iface.__getattribute__(w_name).set_text(phy.__getattribute__(w_name))

    # load plot settings
    [iface.__getattribute__(w_name).set_active(phy.__getattribute__(w_name)) for w_name in
     ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
      'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]

    iface.flipspin.props.adjustment.props.value = phy.flip
    iface.distspin.props.adjustment.props.value = phy.dist

    if data.genes:
        init_sp_genes(iface.sp_genes, data.genes, phy.sel_gene)
    else:
        raise ValueError('no genes')


# from ab12phylo.py

def _exclude(ex_motifs, leaves):
    # filter leaves with exclusion motifs
    if ex_motifs != '':
        # compile exclusion regex.
        queries = ['.*' + re.escape(word.strip()) + '.*' for word in ex_motifs.split(',')]
        if '.*.*' in queries and len(queries) > 1:
            queries.remove('.*.*')
        regex = re.compile(r'|'.join(queries))

        # drop excluded tips
        for dropped_tip in list(filter(regex.match, leaves)):
            leaves.remove(dropped_tip)
    return leaves


def _to_files(leaves):
    # translate to file paths
    with open('metadata.tsv', 'r') as tsv:
        sam_to_file = {leaves.index(sample): _file for sample, _file
                       in [line.split('\t')[0:3:2]  # get first and third column
                           for line in tsv.readlines()] if sample in leaves}
    return [_file for sample, _file in sorted(sam_to_file.items())]


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


def _diversity_stats(gene_now, records, limits, poly):
    # translate dictionary of sequences to numpy int array
    seqs = np.array([seq for seq in records.values()])
    n_sequences = len(records)
    bases = ['A', 'C', 'G', 'T', '-']
    weird_chars = np.setdiff1d(np.unique(seqs), bases)
    codes = dict(zip(bases + list(weird_chars), range(len(bases) + len(weird_chars))))
    coded_seqs = np.vectorize(codes.get)(seqs)

    conserved, biallelic, polyallelic, unknown, gaps = 0, 0, 0, 0, 0
    drop = set()

    # iterate over MSA columns, search for segregating sites
    for j in range(coded_seqs.shape[1]):
        col = coded_seqs[:, j]
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
            if max(col) > 4 and len(col[col > 4]) > limits[1] * n_sequences:
                # there are too many unknown chars in the column
                unknown += 1
                drop.add(j)
            elif max(col) == 4 and len(col[col >= 4]) > limits[0] * n_sequences:
                # there are too many gaps at the site
                gaps += 1
                drop.add(j)
            else:
                # if gaps or unknown characters are at the site, replace them with the most common nucleotide
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
    # print('GENE:%s\tallsites:%d\tconserved:%d\tbiallelic:%d\tpolyallelic:%d\tunknown:%d\tgaps:%d\tpoly:%s\ndrop:%d\n'
    #       % (gene_now, coded_seqs.shape[1], conserved, biallelic, polyallelic, unknown, gaps, poly, len(drop)))

    seg_sites = biallelic
    n_sites = coded_seqs.shape[1] - unknown - gaps
    if poly:
        seg_sites += polyallelic
    else:
        n_sites -= polyallelic

    if seg_sites == 0:
        return '<tr><td>%s</td><td>%d</td><td>0</td><td>--</td><td>--</td><td>--</td><td>--</td>' \
               '<td>--</td><td>%d</td><td>%d</td></tr>' % (gene_now, n_sites, gaps, unknown)

    # crop to allowed sites
    coded_seqs = coded_seqs[:, list(set(range(coded_seqs.shape[1])) - drop)]

    k = 0
    # count pairwise differences
    for combi in itertools.combinations(range(n_sequences), 2):
        k += np.sum(coded_seqs[combi[0]] != coded_seqs[combi[1]])

    # binomial coefficient, k = k^ in formula (10) from Tajima1989
    k = k * 2 / n_sequences / (n_sequences - 1)
    pi = k / n_sites

    a1 = _h(n_sequences)
    theta_w = seg_sites / a1 / n_sites  # per site from Yang2014

    # Tajima's D Formula:
    # d = k - seg_sites/a1 = k - theta_W
    # D = d / sqrt(Var(d))
    try:
        a2 = _qh(n_sequences)
        e1 = 1 / a1 * ((n_sequences + 1) / (3 * n_sequences - 3) - 1 / a1)
        b2 = (2 * (n_sequences * n_sequences + n_sequences + 3)) / (9 * n_sequences * (n_sequences - 1))
        c2 = b2 - ((n_sequences + 2) / (n_sequences * a1)) + (a2 / (a1 * a1))
        e2 = c2 / (a1 * a1 + a2)
        tajima = (k - theta_w) / np.sqrt(e1 * seg_sites + e2 * seg_sites * (seg_sites - 1))
    except (ZeroDivisionError, RuntimeWarning) as z:
        tajima = float('inf')

    n_genotypes = len(np.unique(coded_seqs, axis=0))

    return '<tr><td>%s</td><td>%d</td><td>%d</td><td>%.5f</td><td>%.5f</td>' \
           '<td>%.5f</td><td>%.5f</td><td>%d</td><td>%d</td><td>%d</td></tr>' \
           % (gene_now, n_sites, seg_sites, k, pi, theta_w, tajima, n_genotypes, gaps, unknown)
