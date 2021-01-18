# 2020 Leo Kaindl

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
import toyplot
import toytree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from toyplot import html, pdf, png, svg, color

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo_gui import static, shared
from ab12phylo_gui.static import PATHS as p, seqtoint, SEP, tohex, \
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

    iface.save_plot.connect('clicked', lambda *args: on_save_plot(gui))

    iface.view_msa_ids.set_model(data.tree_anno_model)
    col = Gtk.TreeViewColumn(title='id', cell_renderer=Gtk.CellRendererText(), text=0)
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
    sel = iface.view_msa_ids.get_selection()
    sel.set_mode(Gtk.SelectionMode.MULTIPLE)
    sel.connect('changed', shared.keep_visible,
                iface.parallel_tree.props.vadjustment.props, iface.tempspace)
    iface.msa_eventbox.connect('button_press_event', shared.select_seqs, PAGE, iface.zoomer,
                               iface.view_msa_ids, iface.tempspace)  # in-preview selection
    iface.msa_eventbox.connect('scroll-event', shared.xy_scale, gui, PAGE)  # zooming

    iface.tree = None


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    reload_ui_state(gui)
    shared.load_colorbar(iface.palplot1, gui.wd)
    init_sp_genes(iface.sp_genes, data.genes, phy.sel_gene)
    if shared.get_changed(gui, PAGE):
        start_phy(gui)


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


def start_popgen(widget, gui):
    """Prepare the plotting Thread"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return
    phy.mode = widget.get_name()
    LOG.debug('start_popgen', phy.mode)

    data.pop_model.clear()

    # TODO popgen calc
    init_popgen_columns(iface.popgen)
    return


def init_popgen_columns(tv):
    if len(tv.get_model()) == 0 or tv.get_n_columns > 0:
        return
    for i, ti in enumerate(['gene', 'sites', 'S', 'k', 'π', 'θ',
                            'Tajima\'s D', 'unique', 'gap', 'unknown']):
        tv.append_column(Gtk.TreeViewColumn(
            title=ti, cell_renderer=Gtk.CellRendererText(), text=i))


def on_save_plot(gui):
    """Save the current plot in different formats"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    LOG.debug('rendering tree')
    if phy.png:
        iface.text = 'rendering png'
        png.render(iface.canvas, str(gui.wd / (phy.tx + '.png')), scale=1.6)
    if phy.pdf:
        iface.text = 'rendering pdf'
        pdf.render(iface.canvas, str(gui.wd / (phy.tx + '.pdf')))
    if phy.svg:
        iface.text = 'rendering svg'
        svg.render(iface.canvas, str(gui.wd / (phy.tx + '.svg')))
    if phy.html:
        iface.text = 'rendering html'
        html.render(iface.canvas, str(gui.wd / (phy.tx + '.html')))
    return


def start_phy(gui):
    """Get settings and block GUI"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return
    LOG.debug('start_plotting')

    # parse plot settings only when needed
    [phy.__setattr__(w_name, iface.__getattribute__(w_name).get_active()) for w_name in
     ['rect', 'circ', 'unro', 'tbe', 'fbp', 'supp', 'spec',
      'axis', 'align', 'pmsa', 'pdf', 'svg', 'png', 'nwk', 'html']]
    phy.flip = iface.flipspin.props.adjustment.props.value
    phy.dist = iface.distspin.props.adjustment.props.value
    phy.sel_gene = iface.sp_genes.get_active_text()

    # re-direct to thread
    iface.thread = threading.Thread(target=do_phy1, args=[gui])
    iface.k = (phy.rect + phy.circ + phy.unro) * 2 + 9
    iface.i = 0
    iface.text = 'read tree'
    sleep(.1)
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    return
    # return to main loop


def do_phy1(gui):
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    tree_file = gui.wd / p.fbp if phy.fbp else gui.wd / p.tbe

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
        for t in phy.tips:
            rows = [d.loc[t] for d in dfs]
            if pid:
                phy.ndf[t] = [dict(r) for r in rows if r.pid == max([q.pid for q in rows])][0]
            else:
                phy.ndf[t] = [dict(r) for r in rows][0]

    iface.text = 'fill id column'
    data.tree_anno_model.clear()
    iface.i += 1
    GObject.idle_add(do_phy2, gui)
    sleep(.1)
    return True


def do_phy2(gui):
    """Bounce to idle to size the id column"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.thread.join()
    LOG.info('bouncing phylo thread #1')

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
        c = iface.BLUE if 'accession' in nd else iface.FG
        data.tree_anno_model.append([t, sp, c, iface.BG])

    gui.win.show_all()
    sleep(.1)

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
    phy.g_lens = [(gene, len(seq)) for gene, seq in
                  zip(data.genes, records[next(iter(records))].split(SEP))]
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
    # get model height
    phy.height = iface.view_msa_ids.get_allocated_height() + 50 * phy.axis

    # pre-calculate plot attributes
    # internal node size -> support
    min_size = .2
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

        on_save_plot(gui)
        iface.i += 1

    except ET.ParseError as ex:
        e = 'XML ParseError. Invalid character in a sample ID? Please check metadata.tsv'
        LOG.error(e)
        errors.append(e)
    except ValueError as ex:
        LOG.exception(ex)
        errors.append(ex)

    iface.text = 'placing tree'
    GObject.idle_add(do_phy4, gui, errors)
    sleep(.1)
    return True


def do_phy4(gui, errors):
    """Bounce the plotting to the main loop"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    iface.thread.join()
    LOG.info('bouncing phylo thread #2')
    shared.update(iface, PAGE)
    gui.win.show_all()
    if errors:
        shared.show_notification(gui, 'Errors during tree plotting', errors)
        return
    sleep(.1)

    shared.load_image(iface.zoomer, PAGE, iface.tree_left_vp,
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
    array = np.ma.masked_where(array > static.toint('else'), array)

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
    Path.mkdir(gui.wd / p.phylo.parent, exist_ok=True)
    f.savefig(gui.wd / p.phylo, transparent=True,
              dpi=scale * static.DPI, bbox_inches='tight', pad_inches=0.00001)
    f.savefig(gui.wd / (str(p.phylo)[:-4] + '.svg'), transparent=True,
              dpi=scale * static.DPI, bbox_inches='tight', pad_inches=0.00001)
    plt.close(f)

    shape = array.shape[1], phy.height

    if iface.rasterize.props.active:
        iface.text = 'place PNG MSA: %d:%d' % (array.shape)
        LOG.debug(iface.text)
        sleep(.05)
        shared.load_image(iface.zoomer, PAGE, iface.msa_eventbox, gui.wd / p.phylo,
                          w=shape[0] * 1.6 * shared.get_hadj(iface), h=shape[1] - 50 * phy.axis)
    else:
        iface.text = 'place vector'
        LOG.debug(iface.text)
        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
        canvas.set_size_request(shape[0] * 1.6 * shared.get_hadj(iface), shape[1] - 50 * phy.axis)
        try:
            iface.qal_eventbox.get_child().destroy()
            iface.qal_eventbox.add(canvas)
        except Gtk.Error as ex:
            LOG.error(ex)
    iface.i += 1
    shared.load_colorbar(iface.palplot1, gui.wd)

    for wi in [iface.tree_right, iface.view_msa_ids]:
        wi.set_margin_bottom(50 * phy.axis)

    # save for re-use in popgen part?
    phy.array = array

    iface.text = 'idle'
    sleep(.1)
    GObject.idle_add(stop_phy, gui)
    return True


def stop_phy(gui):
    gui.iface.thread.join()
    LOG.info('phylo thread idle')
    gui.win.show_all()
    gui.data.change_indicator[PAGE] = False
    return False
