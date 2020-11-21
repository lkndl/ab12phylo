# 2020 Leo Kaindl

import gi
from pathlib import Path
import logging
import sys
import hashlib

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from GUI.gtk3 import files, regex, quality, align

LOG = logging.getLogger(__name__)
TOOLS = Path(__file__).resolve().parents[2] / 'ab12phylo' / 'tools'
USER = 'leo.kaindl@tum.de'
SEP = 'SSSSSSSSSS'
RAW_MSA = 'raw_msa.fasta'
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

# TODO continue
# re-fresh page content. automatically called
REFRESH = [files.refresh, regex.refresh, quality.refresh, align.refresh]
# re-run background threads. -> "REFRESH" button
RERUN = {1: regex.start_read, 2: quality.start_trim}
# where the gene selector is visible
SELECT = {2, 4}

# LAMBDAS
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
    gui.iface.reveal_list.props.visible = items
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
    if gui is None:
        gui = widget
    data, iface = gui.data, gui.iface
    page = iface.notebook.get_current_page() if not page else page
    # first integrate changes to the dataset
    if get_changed(gui, page):
        if page == 0:
            regex.reset(gui, do_parse=True)
        elif page == 1:
            if get_errors(gui, page):
                show_message_dialog('There are still errors on the page!')
                return
            if 1 < sum(iface.rx_fired) < 5:
                show_message_dialog('Make sure all columns have been parsed.')
                return
            regex.start_read(gui, run_after=proceed)
            return  # leave this alone
        elif page == 2:
            quality.trim_all(gui)
        elif page == 3:
            print('oh noes')
        set_changed(gui, page, False)

    # then proceed
    iface.notebook.next_page()
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
    RERUN[gui.iface.notebook.get_current_page()](gui)


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
        quality.parse(gui.iface.gene_roll, None, gui)
    # TODO continue


def get_column(list_store, col_idx):
    """Extract a column from a Gtk.ListStore, because they are annoying."""
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col


def update(iface, page):
    """
    Keep the progress bar up-to-date
    :param page: the page to freeze
    """
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
        try:
            rows = [list(row) for row in self]
            coltypes = [type(c) for c in rows[0]]
            return _unpickle_liststore, (self.__class__, coltypes, rows)
        except IndexError as ex:
            cols = self.get_n_columns()
            if cols == 4:
                coltypes = [str, str, str, str]
            elif cols == 3:
                coltypes = [str, bool, bool]
            elif cols == 7:
                coltypes = [str, str, str, str, str, bool, bool, str]
            else:
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
