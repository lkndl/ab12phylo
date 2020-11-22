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


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    iface.gbl = Namespace()

    iface.gaps.connect('changed', lambda *args: iface.gbl.
                       __setattr__('gaps', iface.gaps.get_active_text()))


def refresh(gui):
    """Re-view the page. Load .png images for quick reponsivity."""
    data, iface = gui.data, gui.iface
    if 0 in data.gbl_shape or not (gui.wd / shared.LEFT).exists() or not (gui.wd / shared.RIGHT).exists():
        if not (gui.wd / shared.RAW_MSA).exists():
            return

        # fetch the MSA shape
        i = 0
        for r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta'):
            i += 1
        data.gbl_shape = len(r), i, 0, 0
        iface.msa_shape.set_text('%d : %d' % data.gbl_shape[:2])
        # iface.conserved.set_ad

        return

    # place the png preview
    for wi, image in zip([iface.gbl_left, iface.gbl_right], [gui.wd / shared.LEFT, gui.wd / shared.RIGHT]):
        [child.destroy() for child in wi.get_children()]
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            str(image), width=data.gbl_shape[0], height=data.gbl_shape[1],
            preserve_aspect_ratio=False)
        wi.add(Gtk.Image.new_from_pixbuf(pixbuf))
    # link and resize scrollbar
    iface.gbl_scroll.do_move_slider(iface.gbl_scroll, Gtk.ScrollType.STEP_RIGHT)
    gui.win.show_all()


def prep_gbl(gui):
    pass


def start_gbl(widget, gui):
    """Set-up the Gblocks thread."""
    data, iface = gui.data, gui.iface
    if not data.genes or iface.running or iface.notebook.get_current_page() != PAGE:
        LOG.debug('abort Gblocks')
        return
    iface.frac = 0
    iface.txt = ''
    i = 0
    k = 5


def do_gbl(gui):
    """Run the Gblocks thread."""
    data, iface = gui.data, gui.iface
    iface.txt = 'read pre-trim MSA'
    left_array = np.array([shared.seqtoint(r.seq) for r in SeqIO.parse(gui.wd / shared.RAW_MSA, 'fasta')])

    pass


def stop_gbl(gui):
    """Finish the Gblocks thread"""
    pass
