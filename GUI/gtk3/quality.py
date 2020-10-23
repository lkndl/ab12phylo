import logging
import queue, string
import pandas as pd
from pathlib import Path
from time import sleep
from threading import Thread

import gi
import re
import numpy as np
import seaborn as sns

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GLib

from GUI.gtk3 import commons

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 2
VARIABLE = []


def init(gui):
    data, iface = gui.data, gui.interface
    iface.read_prog.set_visible(False)
    iface.gene_roll.set_model(Gtk.ListStore(str))
    iface.gene_roll.connect('changed', redraw, data, iface)

    iface.accept_rev.set_active(iface.search_rev)
    iface.accept_nophred.set_active(True)

    iface.min_phred.set_adjustment(Gtk.Adjustment(value=30, upper=60, lower=0,
                                                  step_increment=1, page_increment=1))
    iface.min_phred.set_numeric(True)
    iface.min_phred.set_update_policy(Gtk.SpinButtonUpdatePolicy.IF_VALID)

    iface.q_params = dict()

    for w_name in ['min_phred', 'trim_out', 'trim_of', 'bad_stretch']:
        wi = iface.__getattribute__(w_name)
        iface.q_params[w_name] = int(wi.get_text())
        wi.connect('changed', edit, data, iface)
        wi.connect('focus_out_event', parse, data, iface)
        wi.connect('activate', parse, None, data, iface)

    for wi in [iface.trim_out, iface.trim_of, iface.bad_stretch]:
        wi.connect('key-press-event', keypress, data, iface)

    for widget in [iface.accept_rev, iface.accept_nophred]:
        widget.connect('toggled', redraw, data, iface)

    iface.quality_next.connect('clicked', commons.proceed, gui)
    iface.quality_back.connect('clicked', commons.step_back, gui)

    reset(gui)


def reset(gui):
    if len(gui.data.rx_model) > 0:
        data, iface = gui.data, gui.interface
        # queue to share data between threads
        gui.queue = queue.Queue()

        # install timer event to check for new data from the thread
        GLib.timeout_add(50, _on_timer, gui)
        gui.reader = reader_thread(gui)
        gui.reader.start()


class reader_thread(Thread):
    def __init__(self, gui):
        Thread.__init__(self)
        self._gui = gui
        self._rx_model = [row for row in gui.data.rx_model]
        self._wp_model = [row for row in gui.data.wp_model]
        self._csvs = gui.data.csvs
        self._page = gui.notebook.get_children()[PAGE]
        self._queue = gui.queue
        self._show = commons.show_message_dialog
        self._LOG = LOG
        self._csvs = dict()

    def run(self):
        self._page.set_sensitive(False)
        done = 0
        df, box, = None, None
        # read in wellsplates
        LOG.debug('reading wellsplates')
        for row in self._wp_model:
            print(row[:])
            df = pd.read_csv(row[-1], header=None, engine='python')
            df.index = list(range(1, df.shape[0] + 1))
            df.columns = list(string.ascii_uppercase[0:df.shape[1]])
            box = row[0]
            if box in self._csvs:
                self._LOG.error('wellsplate %s already read in. overwrite with %s' % (box, row[-2]))
                self._show(message='wellsplate %s already read in. overwrite with %s' % (box, row[-2]))
            self._csvs[box] = df
            done += 1
            self._queue.put([done, 'plate %s' % row[-2]])

        # read in trace files
        # LOG.debug('reading traces')
        for row in self._rx_model:
            sleep(.3)

            done += 1
            print(done)
            self._queue.put([done, 'trace %s' % row[-2]])

        try:
            self._page.set_sensitive(True)
        except Exception:
            # LOG.error('inactivate notebook page failed')
            pass
        redraw(None, self._gui.data, self._gui.interface)


def _on_timer(gui):
    """updates value of the progress bar"""
    data, iface = gui.data, gui.interface
    tasklen = len(data.rx_model) + len(data.wp_model)
    if not gui.reader.is_alive():
        iface.read_prog.set_visible(False)
        iface.read_prog.set_fraction(0)
        return False

    # # if data available
    # while not gui.queue.empty():
    #     entry = gui.queue.get()
    #     iface.read_prog.set_visible(True)
    #     iface.read_prog.set_fraction(entry[0] / tasklen)
    #     iface.read_prog.set_text('reading %s' % entry[1])

    # keep the timer alive
    return True


def delete_event(widget, event):
    return False


def keypress(widget, event, data, iface):
    key = Gdk.keyval_name(event.keyval)
    if key == 'Up':
        widget.set_text(str(1 + int(widget.get_text())))
        return True
    elif key == 'Down':
        widget.set_text(str(max(0, -1 + int(widget.get_text()))))
        return True


def edit(widget, data, iface):
    LOG.debug('editing ...')
    # filter for numbers only
    value = ''.join([c for c in widget.get_text() if c.isdigit()])
    widget.set_text(value if value else '0')


def parse(widget, event, data, iface):
    LOG.debug('parsing ...')
    # check if the numbers have changed
    pre = iface.q_params[widget.get_name()]
    now = int(widget.get_text())
    delete_event(widget, event)

    # cause redrawing if something changed
    if pre != now:
        iface.q_params[widget.get_name()] = now
        if iface.q_params['trim_out'] > iface.q_params['trim_of']:
            commons.show_message_dialog('cannot draw %d from %d' %
                                        (iface.q_params['trim_out'], iface.q_params['trim_of']))
            widget.set_text('0')
        else:
            redraw(None, data, iface)


def redraw(widget, data, iface, gene='all'):
    LOG.debug('drawing ...')

    # get parameters from interface
    print(iface.q_params)
