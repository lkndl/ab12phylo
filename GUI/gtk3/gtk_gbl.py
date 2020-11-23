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

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, GdkPixbuf

from GUI.gtk3 import shared

LOG = logging.getLogger(__name__)
PAGE = 4

"""The Gblocks page for MSA trimming. Very similar to the sequence trimming page gtk_qal.py"""


# TODO allow deletion at this point!!!!!!


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    iface.view_gbl.set_model(data.gbl_model)
    iface.view_gbl.append_column(Gtk.TreeViewColumn(
        title='id', cell_renderer=Gtk.CellRendererText(), text=0))
    iface.view_gbl.connect('size_allocate', shared.get_dims, iface.gbl_spacer,
                           [iface.gbl_left, iface.gbl_right])

    iface.gbl_preset.set_id_column(0)
    iface.gaps.set_id_column(0)

    iface.gbl = Namespace()
    iface.gbl.params = list()
    for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
        wi = iface.__getattribute__(w_name)
        iface.gbl.__setattr__(w_name, wi.get_adjustment())
        wi.connect('activate', __start_gbl, gui)  # when hitting Enter, re-plot
    for wi in [iface.conserved, iface.flank]:
        wi.connect_after('changed', edit, gui)  # when editing field, check adjustments

    iface.gaps.connect('changed', lambda *args: iface.gbl.  # no assignment in lambda -> use __setattr__
                       __setattr__('gaps', iface.gaps.get_active_text()))
    iface.gbl_preset.connect('changed', re_preset, gui)
    iface.gbl_preset.set_active_id('relaxed')  # trigger initial values. not sure if works?


def refresh(gui):
    """Cause the initial plotting or load .png images for good responsiveness later on."""
    data, iface = gui.data, gui.iface

    if 0 in data.msa_shape or not (gui.wd / shared.LEFT).exists() or not (gui.wd / shared.RIGHT).exists():
        if not (gui.wd / shared.RAW_MSA).exists():
            # stop during initialization
            return

        # start the first re-plotting
        # fetch the MSA shape
        i, r = 0, ''
        data.gbl_model.clear()
        for r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta'):
            data.gbl_model.append([r.id, r.id])
            i += 1
        data.msa_shape = len(r), i, 1, 1  # todo think if this is any use.
        iface.msa_shape.set_text('%d : %d' % data.msa_shape[:2])
        # set the adjustment boundaries of the spin buttons
        re_preset(iface.gbl_preset, gui)
        start_gbl(gui)
        return

    # place the existing png preview
    for wi, image in zip([iface.gbl_left, iface.gbl_right], [gui.wd / shared.LEFT, gui.wd / shared.RIGHT]):
        [child.destroy() for child in wi.get_children()]
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            str(image), width=data.gbl_shape[0], height=data.gbl_shape[1],
            preserve_aspect_ratio=False)
        wi.add(Gtk.Image.new_from_pixbuf(pixbuf))
    gui.win.show_all()


def re_preset(gbl_preset, gui):
    """
    Apply a preset on toggling it.
    :param gbl_preset: the GtkComboBoxText with the presets
    """
    mode = gbl_preset.get_active_text()
    gbl = gui.iface.gbl
    n_sites, n_seqs = gui.data.msa_shape[:2]

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


def __start_gbl(widget, gui):
    start_gbl(gui)


def start_gbl(gui, run_after=None):
    """Set-up the Gblocks thread. This cannot be reached if Gblocks shall be skipped. """
    data, iface = gui.data, gui.iface
    if not data.genes:
        LOG.debug('abort Gblocks')
        return
    elif iface.running:
        shared.show_notification(gui, 'Thread running')
        return
    elif run_after and (gui.wd / shared.MSA).exists() \
            and not shared.get_errors(gui, PAGE):
        shared.set_changed(gui, PAGE, False)
        [do_func(gui) for do_func in run_after]
        return
    elif iface.gbl_preset.get_active_text() == 'skip':
        shared.set_changed(gui, PAGE, False)
        shared.set_errors(gui, PAGE, False)
        shutil.copyfile(shared.RAW_MSA, shared.MSA)
        LOG.info('copied untrimmed MSA')
        return

    # get arguments from adjustments values
    gbl = iface.gbl
    cons, flank, good, bad = [adj.get_value() for adj in
                              [gbl.conserved, gbl.flank, gbl.good_block, gbl.bad_block]]
    gaps = iface.gaps.get_active_text()[0]
    # failsave
    flank = min(cons, flank)
    params = [cons, flank, bad, good, gaps]
    # check if they have changed
    if iface.gbl.params == params:
        shared.show_notification(gui, 'MSA already trimmed, please proceed')
        return
    iface.gbl.params = params

    [[child.destroy() for child in wi.get_children()]
     for wi in [iface.gbl_left, iface.gbl_right]]
    iface.thread = threading.Thread(target=do_gbl, args=[gui])
    iface.run_after = run_after
    iface.running = True
    iface.frac = 0
    iface.txt = ''
    GObject.timeout_add(100, shared.update, iface, PAGE)
    iface.thread.start()
    sleep(.1)
    # return to main loop


def do_gbl(gui):
    """Run the Gblocks thread."""
    data, iface = gui.data, gui.iface

    errors = list()
    iface.frac = .05
    i = 0
    k = 5
    iface.txt = 'read MSA'
    array = np.array([shared.seqtoint(r.seq) for
                      r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta')])
    data.gbl_shape[0] = array.shape[1] * shared.H_SCALER
    # make everything beyond N completely transparent
    array = np.ma.masked_where(array > shared.toint('N'), array)

    i += 1
    iface.frac = i / k
    iface.txt = 'run Gblocks'

    # look for local system Gblocks
    binary = shutil.which('Gblocks')
    local = bool(binary)
    # else pick deployed Gblocks
    binary = binary if binary else shared.TOOLS / 'Gblocks_0.91b' / 'Gblocks'
    LOG.info('%s Gblocks' % ('local' if local else 'packaged'))

    # create base call MARK -t=d sets the mode to nucleotides ... adapt?
    arg = '%s %s -t=d -b1=%d -b2=%d -b3=%d -b4=%d -b5=%s -e=.txt -s=y -p=s; exit 0' \
          % tuple([binary, gui.wd / shared.RAW_MSA] + iface.gbl.params)
    LOG.debug(arg)

    with open(gui.wd / shared.RAW_MSA.parent / 'gblocks.log', 'w') as log_handle:
        try:
            subprocess.run(arg, shell=True, check=True, stdout=log_handle, stderr=log_handle)
        except (OSError, subprocess.CalledProcessError) as e:
            errors.append(str(e))
            log_handle.write(str(e))
            GObject.idle_add(stop_gbl, gui, errors)
            return True

    # parse result
    i += 1
    iface.frac = i / k
    iface.txt = 'parse result'
    shutil.move(gui.wd / shared.RAW_MSA.with_suffix('.fasta.txt'), gui.wd / shared.MSA)
    with open(gui.wd / shared.RAW_MSA.with_suffix('.fasta.txt.txts'), 'r') as fh:
        for line in fh:
            if line.startswith('Flank positions of the') and \
                    line.strip().endswith('selected block(s)'):
                blocks = fh.readline().strip()[9:-1].split(']  [')
                break
    slices = {range(r[0], r[1] + 1) for r in
              [[int(b) for b in block.split('  ')] for block in blocks]}
    LOG.debug(slices)

    i += 1
    iface.frac = i / k
    iface.txt = 'plot MSAs'

    # create a transparency mask
    drop_cols = set()
    while slices:
        drop_cols |= set(slices.pop())
    gbl_mask = np.full(array.shape, shared.ALPHA)
    gbl_mask[:, list(drop_cols)] = 1

    for alpha, win, png_path in zip([1, gbl_mask], [iface.gbl_left, iface.gbl_right],
                                    [shared.LEFT, shared.RIGHT]):
        f = Figure()
        ax = f.add_subplot(111)
        mat = ax.matshow(array, alpha=alpha, cmap=ListedColormap(shared.colors),
                         vmin=-.5, vmax=len(shared.colors) - .5, aspect='auto')
        ax.axis('off')
        f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        f.set_facecolor('none')

        canvas = FigureCanvas(f)
        canvas.set_size_request(data.gbl_shape[0], data.gbl_shape[1])
        try:
            win.add(canvas)
        except Gtk.Error as ex:
            LOG.error(ex)
        pass

        if gui.wd:
            LOG.debug('saving plot')
            Path.mkdir(gui.wd / png_path.parent, exist_ok=True)
            f.savefig(gui.wd / png_path, transparent=True,
                      dpi=600, bbox_inches='tight', pad_inches=0)
        else:
            assert False
        i += 1
        iface.frac = i / k

    # experimentally good point to re-size
    data.gbl_shape[1] = shared.get_dims(iface.view_gbl, None, iface.gbl_spacer,
                                        [iface.gbl_left, iface.gbl_right])
    # re-size
    for wi in [iface.gbl_left, iface.gbl_right]:
        wi.set_max_content_height(data.gbl_shape[1])

    pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(gui.wd / shared.CBAR), 250, 100, preserve_aspect_ratio=True)
    iface.palplot2.set_from_pixbuf(pb)

    sleep(.1)
    GObject.idle_add(stop_gbl, gui, errors)
    return True


def stop_gbl(gui, errors):
    """Finish the Gblocks thread"""
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    gui.win.show_all()
    iface.prog_bar.props.text = 'idle'
    LOG.info('gbl thread idle')
    shared.set_errors(gui, PAGE, bool(errors))
    shared.set_changed(gui, PAGE, False)
    if errors:
        shared.show_notification(gui, 'Errors during MSA trimming', errors)
        return
    if iface.run_after:
        [do_func(gui) for do_func in iface.run_after]
    return
