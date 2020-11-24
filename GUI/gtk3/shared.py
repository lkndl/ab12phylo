# 2020 Leo Kaindl

import gi
from pathlib import Path
import logging
import sys
import copy
import hashlib

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from GUI.gtk3 import gtk_io, gtk_rgx, gtk_qal, gtk_msa, gtk_gbl

LOG = logging.getLogger(__name__)
TOOLS = Path(__file__).resolve().parents[2] / 'ab12phylo' / 'tools'
USER = 'leo.kaindl@tum.de'
SEP = 'SSSSSSSSSS'
RAW_MSA = Path('Trim') / 'raw_msa.fasta'
MSA = 'msa.fasta'
MISSING = 'missing_samples.tsv'
PREVIEW = Path('Trim') / 'trim_preview.png'
CBAR = Path('Trim') / 'colorbar.png'
LEFT = Path('Trim') / 'msa_pretrim.png'
RIGHT = Path('Trim') / 'msa_posttrim.png'
ALPHA = .25
H_SCALER = 3
BUF_SIZE = 128 * 1024
NUCLEOTIDES = ['A', 'C', 'G', 'T', 'N', 'else', '-', ' ', 'S']
KXLIN = {
    'A': (0.92, 1, 0.4, 1),
    'C': (0.46, 1, 0.44, 1),
    'G': (0.16, 0.44, 0.8, 1),
    'T': (1, 0.47, 0.66, 1),
    'N': (0.84, 0.84, 0.84, 0.6),
    'else': (1, 0, 0, 1),
    '-': (1, 1, 1, 0),
    ' ': (1, 1, 1, 0),
    'S': (1, 1, 1, 0)}

# TODO continue project variables
# re-fresh page content. automatically called
REFRESH = [gtk_io.refresh, gtk_rgx.refresh, gtk_qal.refresh, gtk_msa.refresh, gtk_gbl.refresh]
# re-run background threads. -> "REFRESH" button
RERUN = {1: gtk_rgx.start_read, 2: gtk_qal.start_trim, 4: gtk_gbl.start_gbl}
# where the gene selector is visible
SELECT = {2}
algos = {'MAFFT': 'mafft', 'Clustal Omega': 'clustalo', 'MUSCLE': 'muscle', 'T-Coffee': 'tcoffee',
         'RAxML-NG': 'raxml-ng', 'IQ-Tree': 'iqtree', 'FastTree': 'FastTree'}

# LAMBDAS
toalgo = lambda c: algos[c]
toint_map = dict(zip(NUCLEOTIDES, range(len(NUCLEOTIDES))))
toint = lambda c: toint_map.get(c, 5)
seqtoint = lambda s: list(map(toint, s))

togray_map = dict(zip(NUCLEOTIDES, seqtoint('N') * 5 + list(range(5, 9))))
togray = lambda c: togray_map.get(c, 5)
seqtogray = lambda s: list(map(togray, s))

tohex = lambda c: '#' + ''.join([(hex(min(255, int(round(a * 256))))[2:] + '0')[:2].upper() for a in c])
colors = list(map(tohex, map(KXLIN.get, NUCLEOTIDES)))


def get_errors(gui, page):
    return gui.data.errors_indicator[page]


def set_errors(gui, page, errors):
    gui.data.errors_indicator[page] = errors


def get_changed(gui, page):
    return gui.data.change_indicator[page]


def set_changed(gui, page, changed=True):
    """
    Set page to changed and disable further stages; set page to unchanged and enable next page
    :param gui:
    :param page: the index of the current page
    :param changed: toggle
    :return:
    """
    data, iface = gui.data, gui.iface
    if changed:
        # disable later pages
        [page.set_sensitive(False) for page in iface.notebook.get_children()[page + 1: iface.notebook.get_n_pages()]]
        # set later pages to 'changed'
        data.change_indicator[page:] = [True] * (iface.notebook.get_n_pages() - page)
    else:
        # enable the next page
        iface.notebook.get_children()[page + 1].set_sensitive(True)
        # set earlier pages to 'unchanged'
        data.change_indicator[:page + 1] = [False] * (page + 1)


def show_message_dialog(message, list_to_print=None):
    dialog = Gtk.MessageDialog(transient_for=None, flags=0,
                               message_type=Gtk.MessageType.WARNING, text=message)
    # dialog.format_secondary_text('\n'.join(not_found))  # not copyable -> not user-friendly
    if list_to_print:
        txt_buf = Gtk.TextBuffer()
        txt_buf.props.text = '\n'.join(list_to_print)
        txv = Gtk.TextView().new_with_buffer(txt_buf)
        txv.props.wrap_mode = Gtk.WrapMode.WORD_CHAR
        txv.props.margin_end = 20
        txv.props.margin_start = 20
        dialog.get_content_area().add(txv)
    dialog.add_button('OK', 0)
    dialog.get_widget_for_response(0).grab_focus()
    dialog.show_all()  # important
    dialog.run()
    dialog.destroy()
    LOG.debug('showed a message')


def show_notification(gui, msg, items=None):
    gui.iface.reveal_title.set_text(msg)
    gui.iface.reveal_list.props.parent.props.visible = items
    if items:
        gui.iface.reveal_list.props.buffer.props.text = '\n'.join(items)
    gui.iface.revealer.set_reveal_child(True)


def tv_keypress(widget, event, gui, page, selection):
    if Gdk.keyval_name(event.keyval) == 'Delete':
        LOG.debug('delete selected row')
        delete_rows(widget, gui, page, selection)


def delete_rows(widget, gui, page, selection, delete_all=False):
    data, iface = gui.data, gui.iface
    if delete_all:
        delete_all.clear()
    else:
        model, iterator = selection.get_selected_rows()
        [model.remove(model.get_iter(row)) for row in reversed(sorted(iterator))]
    set_changed(gui, page, True)
    REFRESH[page](gui)


def file_hash(file_path):
    h = hashlib.sha256()
    b = bytearray(BUF_SIZE)
    mv = memoryview(b)

    with open(file_path, 'rb', buffering=0) as f:
        for n in iter(lambda: f.readinto(mv), 0):
            h.update(mv[:n])
    return h.hexdigest()


def proceed(widget, gui=None, page=None):
    """
    The function connected to the _Next button. For pages with a background thread,
    this will start it and instruct it to re-run this function afterwards; to make the
    application proceed only upon thread completion.
    # MARK gtk_blast.py will be different and gtk_io.py also is
    :param widget:
    :param gui:
    :param page:
    :return:
    """
    if gui is None:
        gui = widget
    data, iface = gui.data, gui.iface
    page = iface.notebook.get_current_page() if not page else page

    if get_errors(gui, page):
        show_notification(gui, 'There are still errors on the page!')
        return

    if get_changed(gui, page):
        if page == 0:
            gtk_rgx.reset(gui, do_parse=True)
        elif page == 1:
            if 1 < sum(iface.rx_fired) < 5:
                show_notification(gui, 'Make sure all columns have been parsed.')
                return
            gtk_rgx.start_read(gui, run_after=[proceed])
            return  # leave this alone
        elif page == 2:
            gtk_qal.trim_all(gui, run_after=[proceed])
            return
        elif page == 3:
            gtk_msa.start_align(None, gui, run_after=[proceed])
            return
        elif page == 4:
            gtk_gbl.start_gbl(None, gui)
        set_changed(gui, page, False)

    # then proceed
    iface.notebook.next_page()
    # # hide old notifications
    # gui.iface.revealer.set_reveal_child(False)
    LOG.debug('proceeded to page %d' % iface.notebook.get_current_page())


def step_back(widget, gui):
    """Go back to the previous page; set it to un-changed; re-view."""
    data, iface = gui.data, gui.iface
    iface.notebook.prev_page()
    data.page = iface.notebook.get_current_page()
    set_changed(gui, data.page, False)
    REFRESH[data.page](gui)
    LOG.debug('stepped back to page %d' % data.page)


def re_run(gui, *args):
    """Depending on the currently visible page, re-run the matching background task."""
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    # gui.iface.refresh.grab_focus()
    RERUN[page](gui)


def refresh(gui, *args):
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    REFRESH[page](gui)
    # hide or show these two actions depending on applicability
    gui.iface.refresh.props.visible = bool(page in RERUN)
    gui.iface.gene_roll.props.visible = bool(page in SELECT)


def select(gui, *args):
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    if page == 2:  # trim preview
        gtk_qal.parse(gui.iface.gene_roll, None, gui)
    # TODO continue for new pages


def get_column(list_store, col_idx):
    """Extract a column from a Gtk.ListStore, because they are annoying."""
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col


def get_dims(widget, event, spacer, scroll_wins, lower=0):
    """
    Adjust the height of a plot in a GtkScrolledWindow depending on the height of the
    associated labeling column. Adjust the width of the spacer below the labels so that
    the scrollbar below the plot has exactly the width of the latter.
    :param widget: the label column whose re-sizing causes this event
    :param event: will be ignored
    :param spacer: below the label column and next to the scrollbar that will be resized, too
    :param scroll_wins: one or two GtkScrolledWindows containing plots
    """
    w, h = widget.get_allocated_width(), widget.get_allocated_height()
    spacer.set_size_request(w, -1)
    if lower:
        h = min(lower, h)
    for sw in scroll_wins:
        sw.set_max_content_height(h)
        try:
            sw.get_children()[0].set_size_request(w, h)
        except IndexError:
            pass
    return h


def update(iface, page):
    """
    Keep the progress bar up-to-date, and slowly moving rightwards
    :param page: the index of the current page, which will be frozen
    """
    iface.frac = min(iface.frac + 0.001, 1)
    if iface.running:
        iface.notebook.get_children()[page].set_sensitive(False)
        iface.prog_bar.set_visible(True)
        iface.prog_bar.set_fraction(iface.frac)
        iface.prog_bar.set_text(iface.txt)
        return True
    else:
        iface.notebook.get_children()[page].set_sensitive(True)
        iface.prog_bar.set_visible(False)
        return False


def bind_accelerator(accelerators, widget, accelerator, signal='clicked'):
    key, mod = Gtk.accelerator_parse(accelerator)
    widget.add_accelerator(signal, accelerators, key, mod, Gtk.AccelFlags.VISIBLE)


class picklable_liststore(Gtk.ListStore):
    """
    kudos go to samplebias on https://stackoverflow.com/a/5969700
    """

    def __reduce__(self):
        rows = list()
        try:
            rows = [list(row) for row in self]
            coltypes = [type(c) for c in rows[0]]
            return _unpickle_liststore, (self.__class__, coltypes, rows)
        except IndexError as ex:
            cols = self.get_n_columns()
            # allow saving of emptry data stores
            if cols == 2:
                coltypes = [str, str]
            elif cols == 4:
                coltypes = [str, str, str, str]
            elif cols == 3:
                coltypes = [str, bool, bool]
            elif cols == 7:
                coltypes = [str, str, str, str, str, bool, bool, str]
            else:
                print('assert False')
                assert False
            return _unpickle_liststore, (self.__class__, coltypes, rows)
        except Exception as ex:
            LOG.exception(ex)


def _unpickle_liststore(cls, col_types, rows):
    inst = cls.__new__(cls)
    inst.__init__(*col_types)
    for row in rows:
        inst.append(row)
    return inst