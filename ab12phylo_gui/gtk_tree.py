# 2020 Leo Kaindl

import json
import logging
import threading
from pathlib import Path
from time import sleep

import gi
from Bio import SeqIO


gi.require_version('Gtk', '3.0')

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 7


def init(gui):
    """Initialize the page. Connect buttons"""
    data, iface = gui.data, gui.iface
    pass


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    data, iface = gui.data, gui.iface
    # algo = static.toalgo(iface.ml_stack.get_visible_child_name())
    # check_MSA(None, gui, algo)
    pass

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
