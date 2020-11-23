# 2020 Leo Kaindl


import logging
import shutil
import subprocess
import threading
import numpy as np
from Bio import SeqIO
from argparse import Namespace

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, GdkPixbuf

from GUI.gtk3 import shared
from ab12phylo import msa

LOG = logging.getLogger(__name__)
PAGE = 4

"""The Gblocks page for MSA trimming. Very similar to the sequence trimming page gtk_qal.py"""


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    iface.view_gbl.set_model(data.qal_model)
    iface.view_gbl.append_column(Gtk.TreeViewColumn(
        title='id', cell_renderer=Gtk.CellRendererText(), text=0))
    iface.view_gbl.connect('size_allocate', shared.get_dims, iface.gbl_spacer,
                           [iface.gbl_left, iface.gbl_right])

    iface.gbl_preset.set_id_column(0)
    iface.gaps.set_id_column(0)

    iface.gbl = Namespace()
    for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
        wi = iface.__getattribute__(w_name)
        iface.gbl.__setattr__(w_name, wi.get_adjustment())
        wi.connect('activate', start_gbl, gui)  # when hitting Enter, re-plot
    for wi in [iface.conserved, iface.flank]:
        wi.connect('changed', edit, gui)  # when editing field, check adjustments

    iface.gaps.connect('changed', lambda *args: iface.gbl.  # no assignment in lambda -> use __setattr__
                       __setattr__('gaps', iface.gaps.get_active_text()))
    iface.gbl_preset.connect('changed', re_preset, gui)

    iface.gbl_preset.set_active_id('relaxed')  # trigger initial values


def refresh(gui):
    """Re-view the page. Load .png images for good responsiveness."""
    data, iface = gui.data, gui.iface
    if 0 in data.msa_shape or not (gui.wd / shared.LEFT).exists() or not (gui.wd / shared.RIGHT).exists():
        if not (gui.wd / shared.RAW_MSA).exists():
            # stop during initialization
            return

        # fetch the MSA shape
        i, r = 0, ''
        for r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta'):
            i += 1
        data.msa_shape = len(r), i, 0, 0
        iface.msa_shape.set_text('%d : %d' % data.msa_shape[:2])
        # set the adjustment boundaries of the spin buttons
        re_preset(iface.gbl_preset, gui)
        return

    # place the png preview
    data.gbl_shape[1] = shared.get_dims(iface.view_gbl, None, iface.gbl_spacer,
                                        [iface.gbl_left, iface.gbl_right])
    for wi, image in zip([iface.gbl_left, iface.gbl_right], [gui.wd / shared.LEFT, gui.wd / shared.RIGHT]):
        [child.destroy() for child in wi.get_children()]
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            str(image), width=data.gbl_shape[0], height=data.gbl_shape[1],
            preserve_aspect_ratio=False)
        wi.add(Gtk.Image.new_from_pixbuf(pixbuf))
    # link and resize scrollbar
    iface.gbl_scroll.do_move_slider(iface.gbl_scroll, Gtk.ScrollType.STEP_RIGHT)
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
    flank.configure(max(adj.get_value(), flank.get_value()), adj.get_value(), flank.get_upper(), 1, 0, 0)


def start_gbl(widget, gui, run_after=None):
    """Set-up the Gblocks thread."""
    data, iface = gui.data, gui.iface
    if not data.genes or iface.running or iface.notebook.get_current_page() != PAGE:
        LOG.debug('abort Gblocks')
        return
    iface.frac = 0
    iface.txt = ''
    iface.run_after =run_after
    i = 0
    k = 5


def do_gbl(gui):
    """Run the Gblocks thread."""
    data, iface = gui.data, gui.iface
    iface.txt = 'read pre-trim MSA'
    left_array = np.array([shared.seqtoint(r.seq) for r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta')])
    data.gbl_shape[0] = left_array.shape[1] * 4

    data.gbl_shape[1] = shared.get_dims(iface.view_gbl, None, iface.gbl_spacer,
                                        [iface.gbl_left, iface.gbl_right])
    # canvas.set_size...
    pass


def stop_gbl(gui):
    """Finish the Gblocks thread"""
    data, iface = gui.data, gui.iface

    shared.set_changed(gui, PAGE, False)
    if iface.run_after:
        [do_func(gui)for do_func in iface.run_after]
    return
