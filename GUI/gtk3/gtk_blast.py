# 2020 Leo Kaindl

import json
import logging
import re
import string
import sys
import threading
import webbrowser
from pathlib import Path
from time import sleep

import gi
import pandas as pd
import requests, random
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import shared, gtk_qal
from ab12phylo import raxml

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 5




def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface
    pass


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree."""
    pass


def start_BLAST(widget, gui):
    """Set-up the BLAST thread."""
    pass


def do_BLAST(gui):
    """Run BLAST thread."""
    pass


def stop_BLAST(gui):
    """Finish the Gblocks thread"""
    pass
