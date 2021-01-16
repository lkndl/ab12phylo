# 2020 Leo Kaindl

import json
import logging
import threading
import toytree
import toyplot
import numpy as np
from pathlib import Path
from time import sleep

import gi
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from ab12phylo_gui import static, shared

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 7


def init(gui):
    """Initialize the page. Connect buttons"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface

    for w_name in ['gap_share', 'unk_share']:
        iface.__getattribute__(w_name).get_adjustment().connect(
            'value-changed', lambda adj: phy.__setattr__(w_name, adj.props.value))

    for w_name in ['matches', 'subtree']:
        iface.__getattribute__(w_name).connect('clicked', start_pop, gui)


    # iface.view_msa_ids.set_model(phy.phy_model)


def start_pop(widget, gui):
    """Prepare the plotting Thread"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return

    phy.mode = widget.get_name()
    phy.query = iface.query.get_text()
    phy.exclude = iface.exclude.get_text()

    if iface.popgen.get_n_columns() < 10:
        iface.popgen.set_model(data.pop_model)
        for i, ti in enumerate(['gene', 'sites', 'S', 'k', 'π', 'θ',
                                'Tajima\'s D', 'unique', 'gap', 'unknown']):
            iface.popgen.append_column(Gtk.TreeViewColumn(
                title=ti, cell_renderer=Gtk.CellRendererText(), text=i))

    data.pop_model.clear()
    # TODO
    return


def start_phy(gui):
    """Get settings and block GUI"""
    data, phy, iface = gui.data, gui.data.phy, gui.iface
    if iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
        return

    phy.args = [wi.get_name() for wi in
                [iface.rect, iface.circ, iface.unro,
                 iface.to_pdf, iface.to_svg, iface.to_png, iface.to_nwk,
                 iface.supp, iface.spec] if wi.get_active()]



def do_phy(gui):
    """Plot the tree"""
    pass


def do_phyMSA(gui):
    """Plot the MSA"""


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    data, iface = gui.data, gui.iface
    # algo = static.toalgo(iface.ml_stack.get_visible_child_name())
    # check_MSA(None, gui, algo)
    shared.load_colorbar(iface.palplot1, gui.wd)

    #
    # def check_MSA(widget, gui, algo):
    #     pass
    #
    #
    # def start_ML(gui):
    #     """Set-up the Gblocks thread"""
    #     pass
    #
    #
    # def do_ML(gui):
    #     """Run the Gblocks thread"""
    #     pass
    #
    #
    # def stop_ML(gui, errors):
    #     """Finish the ML inference thread"""
    #     iface = gui.iface
    #     iface.running = False
    #     iface.thread.join()
    #     # TODO sensitivize a container?
    #     gui.win.show_all()
    #     iface.prog_bar.props.text = 'idle'
    #     LOG.info('ml thread idle')
    #     shared.set_errors(gui, PAGE, bool(errors))
    #     shared.set_changed(gui, PAGE, False)
    #     if errors:
    #         shared.show_notification(gui, 'Errors during ML inference', errors)
    #         return
    #     if iface.run_after:
    #         [do_func(gui) for do_func in iface.run_after]
    #     return
