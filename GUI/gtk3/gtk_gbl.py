# 2020 Leo Kaindl

import logging
import shutil
import subprocess
import threading
from argparse import Namespace
from pathlib import Path
from time import sleep

import gi
import numpy as np
from Bio import SeqIO
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

import static

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import shared
from static import PATHS

LOG = logging.getLogger(__name__)
PAGE = 4

"""The Gblocks page for MSA trimming. Very similar to the sequence trimming page gtk_qal.py"""


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    iface.view_gbl.set_model(data.gbl_model)
    iface.view_gbl.append_column(Gtk.TreeViewColumn(
        title='id', cell_renderer=Gtk.CellRendererText(), text=0))
    iface.view_gbl.connect('check-resize', shared.get_height_resize, iface.gbl_spacer,
                           [iface.gbl_left, iface.gbl_right])

    iface.gbl_preset.set_id_column(0)
    iface.gaps.set_id_column(0)

    iface.gbl = Namespace()
    iface.gbl.params = list()
    for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
        wi = iface.__getattribute__(w_name)
        iface.gbl.__setattr__(w_name, wi.get_adjustment())
        wi.connect('activate', lambda *args: start_gbl(gui))  # when hitting Enter, re-plot. overshadowed by Refresh
    for wi in [iface.conserved, iface.flank]:
        wi.connect_after('changed', edit, gui)  # when editing field, check adjustments

    iface.gaps.connect('changed', lambda *args: iface.gbl.  # no assignment in lambda -> use __setattr__
                       __setattr__('gaps', iface.gaps.get_active_text()))
    iface.gbl_preset.connect_after('changed', re_preset, gui)
    iface.gbl_preset.set_active_id('balanced')

    sel = iface.view_gbl.get_selection()
    sel.set_mode(Gtk.SelectionMode.MULTIPLE)
    sel.connect('changed', shared.keep_visible,
                iface.parallel_gbl.props.vadjustment.props, iface.gbl)
    iface.view_gbl.connect('key_press_event', shared.delete_and_ignore_rows,
                           gui, PAGE, sel, iface.gbl)  # in-preview deletion
    iface.gbl_eventbox.connect('button_press_event', shared.select_seqs, PAGE, iface.zoom,
                               iface.view_gbl, iface.gbl)  # in-preview selection
    iface.gbl_eventbox.connect('scroll-event', shared.xy_scale, gui, PAGE)  # zooming


def refresh(gui):
    """Cause the initial plotting or load .png images for good responsiveness later on."""
    data, iface = gui.data, gui.iface

    # start the first re-plotting
    if 0 in data.msa_shape or not (gui.wd / PATHS.left).exists() \
            or not (gui.wd / PATHS.right).exists():
        # fetch the MSA shape
        data.gbl_model.clear()
        shared_ids = sorted(set.intersection(*data.gene_ids.values()))
        [data.gbl_model.append([_id]) for _id in shared_ids]
        i = len(shared_ids)
        data.msa_shape = [sum(data.msa_lens), i, -1, i]  # width-height-width_after-height
        iface.msa_shape.set_text('%d : %d' % tuple(data.msa_shape[:2]))
        iface.msa_shape_trimmed.set_text('-1 : %d' % i)
        # set the adjustment boundaries of the spin buttons
        re_preset(iface.gbl_preset, gui)
        start_gbl(gui)
        return

    # else load the pre-existing PNGs
    try:
        iface.msa_shape.set_text('%d : %d' % tuple(data.msa_shape[:2]))
        iface.msa_shape_trimmed.set_text(
            '%d : %d' % (data.msa_shape[2] - (len(data.genes) - 1) * len(static.SEP),
                         data.msa_shape[3]))
    except TypeError:
        pass
    # place the existing png
    x_ratio = data.msa_shape[2] / data.msa_shape[0]
    shared.load_image(iface.zoom, PAGE, iface.gbl_left_vp, gui.wd / PATHS.left,
                      width=data.gbl_shape[0] * shared.get_hadj(iface),
                      height=data.gbl_shape[1])
    shared.load_image(iface.zoom, PAGE, iface.gbl_right_vp, gui.wd / PATHS.right,
                      width=data.gbl_shape[0] * shared.get_hadj(iface) * x_ratio,
                      height=data.gbl_shape[1])
    shared.load_colorbar(iface.palplot2, gui.wd)
    re_preset(iface.gbl_preset, gui)
    gui.win.show_all()


def re_preset(gbl_preset, gui):
    """
    Apply a preset on toggling it.
    :param gbl_preset: the GtkComboBoxText with the presets
    """
    mode = gbl_preset.get_active_text()
    gbl = gui.iface.gbl
    n_sites, n_seqs = gui.data.msa_shape[:2]

    assert type(n_sites) == int
    # block / un-block SpinButtons
    [gui.iface.__getattribute__(w_name).set_sensitive(mode != 'skip')
     for w_name in ['conserved', 'flank', 'good_block', 'bad_block', 'gaps']]

    if mode == 'skip':
        return
    elif mode == 'relaxed':
        conserved = n_seqs // 2 + 1
        flank = conserved
        gaps = 'half'
        good_block = 5
        bad_block = 8
    elif mode == 'balanced':
        conserved = n_seqs // 2 + 1
        flank = min(n_seqs // 4 * 3 + 1, n_seqs)
        gaps = 'half'
        good_block = 5
        bad_block = 8
    elif mode == 'default':
        conserved = n_seqs // 2 + 1
        flank = min(int(n_seqs * 0.85) + 1, n_seqs)
        gaps = 'none'
        good_block = 10
        bad_block = 8
    elif mode == 'strict':
        conserved = int(n_seqs * .9)
        flank = conserved
        gaps = 'none'
        good_block = 10
        bad_block = 8
    else:
        LOG.error('illegal mode: ' + str(mode))
        assert False

    # configure(value, lower, upper, step-increment=1, page-increment=0, page-size=0)
    gbl.conserved.configure(conserved, n_seqs // 2 + 1, n_seqs, 1, 0, 0)
    gbl.flank.configure(flank, n_seqs // 2 + 1, n_seqs, 1, 0, 0)
    gbl.good_block.configure(good_block, 2, n_sites, 1, 0, 0)
    gbl.bad_block.configure(bad_block, 0, n_sites, 1, 0, 0)
    gui.iface.gaps.set_active_id(gaps)


def edit(wi, gui):
    """
    Just make sure flank is never smaller than conserved
    """
    LOG.debug('editing')
    adj = wi.get_adjustment()
    flank = gui.iface.gbl.flank
    conserved = gui.iface.gbl.conserved
    if adj == flank:
        return  # break some loops
    # conserved has changed
    flank.configure(max(adj.get_value(), flank.get_value()),  # value
                    adj.get_value(), flank.get_upper(), 1, 0, 0)  # the others


def start_gbl(gui, run_after=None):
    """Set-up the Gblocks thread. This cannot be reached if Gblocks shall be skipped. """
    data, iface = gui.data, gui.iface
    if not data.genes:
        LOG.debug('abort Gblocks')
        return
    elif iface.running:
        shared.show_notification(gui, 'Thread running')
        return
    elif run_after and (gui.wd / PATHS.msa).exists() \
            and not shared.get_errors(gui, PAGE):
        shared.set_changed(gui, PAGE, False)
        [do_func(gui) for do_func in run_after]
        return

    data.gbl_model.clear()
    if iface.gbl_preset.get_active_text() != 'skip':
        # get arguments from adjustments values
        gbl = iface.gbl  # the param Namespace
        cons, flank, good, bad = [adj.get_value() for adj in
                                  [gbl.conserved, gbl.flank, gbl.good_block, gbl.bad_block]]
        gaps = iface.gaps.get_active_text()[0]
        # failsafe
        if flank < cons:
            flank = cons
        params = [cons, flank, bad, good, gaps]
        # check if they have changed
        if iface.gbl.params == params:
            shared.show_notification(gui, 'MSA already trimmed, please proceed')
            return
        iface.gbl.params = params

    iface.thread = threading.Thread(target=do_gbl, args=[gui])
    iface.run_after = run_after
    iface.running = True
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    sleep(.1)
    # return to main loop


def do_gbl(gui):
    """Run the Gblocks thread."""
    data, iface = gui.data, gui.iface

    errors = list()
    iface.frac = .05
    iface.i = 0
    iface.k = len(data.genes) + 6

    if iface.gbl_preset.get_active_text() != 'skip':
        iface.text = 'prep'
        LOG.debug(iface.text)
        # look for local system Gblocks
        binary = shutil.which('Gblocks')
        local = bool(binary)
        # else pick deployed Gblocks
        binary = binary if binary else static.TOOLS / 'Gblocks_0.91b' / 'Gblocks'
        LOG.info('%s Gblocks' % ('local' if local else 'packaged'))

        # create base call MARK -t=d sets the mode to nucleotides ... adapt?
        arg = '%s %s -t=d -b1=%d -b2=%d -b3=%d -b4=%d -b5=%s -e=.txt -s=y -p=s; exit 0' \
              % tuple([binary, '%s'] + iface.gbl.params)
        LOG.debug(arg)

    # fetch IDs
    shared_ids = sorted(set.intersection(*data.gene_ids.values()))
    array = np.empty(shape=(len(shared_ids), 0), dtype=int)
    [data.gbl_model.append([_id]) for _id in shared_ids]
    msa_lens = list()
    slices = set()
    iface.i += 1

    for gene in data.genes:
        iface.text = '%s: read MSA' % gene
        LOG.debug(iface.text)
        raw_msa = gui.wd / gene / ('%s_raw_msa.fasta' % gene)
        msa = gui.wd / gene / ('%s_msa.fasta' % gene)
        records = {r.id: r.seq.upper() for r in SeqIO.parse(raw_msa, 'fasta')}
        ar = np.array([static.seqtoint(records[_id]) for _id in shared_ids])
        assert ar.shape[0] == len(shared_ids)
        msa_lens.append(ar.shape[1])
        array = np.hstack((array, ar,))  # shared.SEP would need to be stacked here
        data.msa_shape[:2] = array.shape[::-1]

        if iface.gbl_preset.get_active_text() == 'skip':
            shutil.copy(raw_msa, msa)
            iface.i += 1
            continue

        iface.text = '%s: run Gblocks' % gene
        LOG.debug(iface.text)
        with open(gui.wd / gene / 'gblocks.log', 'w') as log_handle:
            try:
                subprocess.run(arg % raw_msa, shell=True, check=True, stdout=log_handle, stderr=log_handle)
            except (OSError, subprocess.CalledProcessError) as e:
                errors.append(str(e))
                log_handle.write(str(e))
                continue

        # parse result
        iface.text = '%s: parse result' % gene
        LOG.debug(iface.text)

        shutil.move(raw_msa.with_suffix('.fasta.txt'), msa)
        with open(raw_msa.with_suffix('.fasta.txt.txts'), 'r') as fh:
            for line in fh:
                if line.startswith('Flank positions of the') and \
                        line.strip().endswith('selected block(s)'):
                    lane = fh.readline()
                    line_blocks = lane.strip()[9:-1].split(']  [')
                    break
        LOG.debug(line_blocks)
        if line_blocks == ['']:
            err = '%s: no good blocks' % gene
            LOG.error(err)
            errors.append(err)
            continue

        shift = sum(msa_lens[:-1])
        slices |= {range(r[0] - 1, r[1]) for r in
                   [[int(b) + shift for b in block.split('  ')] for block in line_blocks]}
        iface.i += 1

    iface.text = 'concatenating MSAs'
    if 'aligner' not in iface:
        iface.aligner, cmd = shared.get_msa_build_cmd(
            static.toalgo(iface.msa_algo.get_active_text()), gui.wd, data.genes)
    iface.aligner.reset_paths(gui.wd, gui.wd / PATHS.msa)
    data.msa_shape[2], data.msa_shape[3] = iface.aligner.concat_msa(gui=shared_ids)
    iface.text = 'computing SHA256 hash'
    LOG.debug(iface.text)
    shared.get_hashes(gui, PATHS.msa, PAGE)

    iface.i += 1
    iface.text = 'plot MSAs'
    LOG.debug(iface.text)

    data.gbl_shape[0] = array.shape[1] * shared.get_hadj(iface)
    # make gaps transparent
    array = np.ma.masked_where(array > static.toint('else'), array)

    # create a transparency mask
    blocks = set()
    while slices:
        blocks |= set(slices.pop())
    blocks = sorted(list(blocks))
    gbl_mask = np.full(array.shape, static.ALPHA)
    gbl_mask[:, blocks] = 1

    data.msa_shape[2] = len(blocks)  # for completeness
    LOG.debug('msa shape: %s' % str(data.msa_shape))
    x_ratio = data.msa_shape[2] / data.msa_shape[0]
    LOG.debug('x ratio: %.3f' % x_ratio)

    # adjust maximum size
    scale = 6
    while max(data.msa_shape[:2]) * scale > 2 ** 14:
        scale -= 1
    LOG.debug('scaling gbl with %d' % scale)

    if iface.gbl_preset.get_active_text() != 'skip':
        for alpha, blocks, gtk_bin, png_path, x_ratio, width \
                in zip([gbl_mask, 1], [range(array.shape[1]), blocks],
                       [iface.gbl_left_vp, iface.gbl_right_vp],
                       [PATHS.left, PATHS.right], [1, x_ratio],
                       [data.msa_shape[0], data.msa_shape[2]]):
            f = Figure()  # figsize=(width / shared.DPI, data.msa_shape[1] / shared.DPI), dpi=shared.DPI)  # figaspect(data.msa_shape[1] / width))
            f.set_facecolor('none')
            f.set_figheight(data.msa_shape[1] / static.DPI)
            f.set_figwidth(max(1, width) / static.DPI)
            ax = f.add_subplot(111)
            ax.axis('off')
            f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
            mat = ax.matshow(array[:, blocks], alpha=alpha, cmap=ListedColormap(static.colors),
                             vmin=-.5, vmax=len(static.colors) - .5, aspect='auto')

            iface.i += 1
            iface.text = 'save PNG'
            LOG.debug(iface.text)
            Path.mkdir(gui.wd / png_path.parent, exist_ok=True)
            f.savefig(gui.wd / png_path, transparent=True,
                      dpi=scale * static.DPI, bbox_inches='tight', pad_inches=0.00001)

            # experimentally good point to re-size
            data.gbl_shape[1] = shared.get_height_resize(iface.view_gbl, None, iface.gbl_spacer,
                                                         [iface.gbl_left, iface.gbl_right])

            if iface.rasterize.props.active:
                iface.text = 'place PNG'
                LOG.debug(iface.text)
                shared.load_image(iface.zoom, PAGE, gtk_bin, gui.wd / png_path,
                                  data.gbl_shape[0] * x_ratio * shared.get_hadj(iface), data.gbl_shape[1])
            else:
                iface.text = 'place vector'
                LOG.debug(iface.text)
                canvas = FigureCanvas(f)
                canvas.set_size_request(len(blocks) * shared.get_hadj(iface), data.gbl_shape[1])  # width, height
                try:
                    gtk_bin.remove(gtk_bin.get_child())
                    gtk_bin.add(canvas)
                except Gtk.Error as ex:
                    LOG.error(ex)
            iface.i += 1

    # re-size
    for wi in [iface.gbl_left, iface.gbl_right]:
        wi.set_max_content_height(data.gbl_shape[1])

    shared.load_colorbar(iface.palplot2, gui.wd)

    iface.text = 'idle'
    iface.frac = 1
    sleep(.1)
    GObject.idle_add(stop_gbl, gui, errors)
    return True


def stop_gbl(gui, errors):
    """Finish the Gblocks thread"""
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    gui.win.show_all()
    LOG.info('gbl thread idle')
    shared.set_errors(gui, PAGE, bool(errors))
    shared.set_changed(gui, PAGE, False)
    iface.msa_shape.set_text('%d : %d' % tuple(gui.data.msa_shape[:2]))
    iface.msa_shape_trimmed.set_text(
        '%d : %d' % (gui.data.msa_shape[2], gui.data.msa_shape[3]))
    if iface.gbl.flank.props.upper != gui.data.msa_shape[1]:
        re_preset(iface.gbl_preset, gui)
    if errors:
        shared.show_notification(gui, 'Errors during MSA trimming', errors)
        return
    if iface.run_after:
        [do_func(gui) for do_func in iface.run_after]
    return


def drop_seqs(gui):
    """
    Overwrite 'old' trimmed MSA file that still contained dropped records.
    :param gui:
    :return:
    """
    if 'ignore_set' not in gui.iface.gbl:
        return
    with open(gui.wd / (PATHS.msa + '_TEMP'), 'w') as fasta:
        for record in SeqIO.parse(gui.wd / PATHS.msa, 'fasta'):
            if record.id not in gui.iface.gbl.ignore_set:
                SeqIO.write(record, fasta, 'fasta')
    LOG.debug('dropped %d sequences from trimmed MSA' % len(gui.iface.gbl.ignore_set))
    del gui.iface.gbl.ignore_set
    shutil.move(gui.wd / (PATHS.msa + '_TEMP'), gui.wd / PATHS.msa)
