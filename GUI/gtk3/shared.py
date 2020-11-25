# 2020 Leo Kaindl

import gi
from pathlib import Path
from time import sleep
import logging
import sys
import copy
import hashlib

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf

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
            gtk_gbl.start_gbl(gui)
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


def load_image(gtk_bin, img_path, width, height):
    LOG.debug('load image')
    child = gtk_bin.get_child()
    if type(child) != Gtk.Image:
        gtk_bin.remove(child)
        child = Gtk.Image()
        gtk_bin.add(child)

    pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(img_path), width=width, height=height,
        preserve_aspect_ratio=False)
    child.set_from_pixbuf(pb)


def load_colorbar(gtk_image, wd):
    gtk_image.set_from_pixbuf(GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(wd / CBAR), width=250, height=100, preserve_aspect_ratio=True))


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
    LOG.debug('re-sizing')
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
    print('new height: %d' % h, file=sys.stderr)
    return h


def hadj(iface):
    return iface.hadj.adj.get_value() * 2


def scale(gtk_bin, factor):
    """Horizontally stretch a pixbuf."""
    child = gtk_bin.get_child()
    if type(child) != Gtk.Image:
        LOG.info('can not rescale %s' % str(type(child)))
        return False
    pb = child.get_pixbuf()
    child.set_from_pixbuf(pb.scale_simple(
        min(14000, pb.get_width() * factor), pb.get_height(),
        GdkPixbuf.InterpType.NEAREST))
    return True


def h_adjust(adj, gui):
    """
    Horizontally scale a preview.
    """
    data, iface = gui.data, gui.iface
    # with adj.hander_block(iface.hadj.handle):
    val = max(.2, adj.get_value())
    page = gui.iface.notebook.get_current_page()
    shift = val - iface.hadj.bak
    if shift > 0 > iface.hadj.trend:
        iface.hadj.trend = shift

    factor = val / iface.hadj.bak
    iface.hadj.trend += factor - 1

    if factor >= 2:
        LOG.debug('re-loading images')
        # re-load and re-set
        iface.hadj.trend = 0
        if page == 2:
            load_image(iface.qal_eventbox, gui.wd / PREVIEW,
                       data.qal_shape[0] * val * 2, data.qal_shape[1])
        elif page == 4:
            load_image(iface.gbl_left_eventbox, gui.wd / LEFT,
                       data.msa_shape[0] * val * 2, data.gbl_shape[1])
            sleep(.1)
            load_image(iface.gbl_right_eventbox, gui.wd / RIGHT,
                       data.msa_shape[2] * val * 2, data.gbl_shape[1])
    else:
        # scale in-place
        LOG.debug('scale: %.2f fold' % factor)
        if page == 2:
            scale(iface.qal_eventbox, factor)
        elif page == 4:
            scale(iface.gbl_left_eventbox, factor)
            sleep(.1)
            scale(iface.gbl_right_eventbox, factor)
    iface.hadj.bak = val
    sleep(.1)


def update(iface, page):
    """
    Keep the progress bar up-to-date, and slowly moving rightwards
    :param page: the index of the current page, which will be frozen
    """
    iface.frac = min(
        max(iface.i / iface.k,  # base progress from caller thread iteration
            iface.frac + 0.001),  # pretend to proceed
        (iface.i + 1) / iface.k,  # but do not pass next iteration level
        1)  # and never pass 1
    if iface.running:
        iface.notebook.get_children()[page].set_sensitive(False)
        iface.prog_bar.set_fraction(iface.frac)
        for wi in [iface.prog_bar, iface.prog_label]:
            wi.set_visible(True)
            wi.set_text(iface.text)
        return True
    else:
        iface.notebook.get_children()[page].set_sensitive(True)
        [wi.set_visible(False) for wi in [iface.prog_bar, iface.prog_label]]
        return False


def bind_accelerator(accelerators, widget, accelerator, signal='clicked'):
    key, mod = Gtk.accelerator_parse(accelerator)
    widget.add_accelerator(signal, accelerators, key, mod, Gtk.AccelFlags.VISIBLE)


def select_seqs(event_box, loc, tv, ns):
    """
    Select sequences from a trim preview directly in the image, in an expected way.
    Shift focus to the labeling treeview on the left, which is already <Delete>-sensitive
    :param event_box:
    :param loc:
    :param sel:
    :param ns:
    :param args:
    :return:
    """
    tv.grab_focus()
    sel = tv.get_selection()

    accel_mask = Gtk.accelerator_get_default_mod_mask()
    rect, baseline = event_box.get_allocated_size()
    mo, tree_path_iterator = sel.get_selected_rows()
    idcs = {tp[0] for tp in tree_path_iterator}
    idx = int(loc.y / rect.height * len(mo))
    if idx == '' or idx < 0:
        sel.unselect_all()
        return
    if loc.state & accel_mask == Gdk.ModifierType.CONTROL_MASK:
        sel.select_path(idx) if idx not in idcs else sel.unselect_path(Gtk.TreePath(idx))
    elif loc.state & accel_mask == Gdk.ModifierType.SHIFT_MASK and 'previous' in ns:
        m, n = min(idx, ns.previous), max(idx, ns.previous)
        new_idcs = set(range(m, n))
        if all(idc in idcs for idc in new_idcs):
            sel.unselect_range(Gtk.TreePath(m), Gtk.TreePath(n))
        else:
            sel.select_range(Gtk.TreePath(m), Gtk.TreePath(n))
    else:
        # select only the one clicked on
        sel.unselect_all()
        sel.select_path(idx)
    ns.previous = idx


def delete_and_ignore_rows(widget, event, gui, page, sel, ns):
    """
    Keep track of the rows that will not be written to the next fasta and delete them from the treeview.
    """
    if Gdk.keyval_name(event.keyval) == 'Delete':
        model, tree_path_iterator = sel.get_selected_rows()
        if 'ignore_set' not in ns:
            ns.ignore_set = set()
        for row in reversed(sorted(tree_path_iterator)):
            ns.ignore_set.add(model[row[:]][0])
            if page == 2:
                model.remove(model.get_iter(row))
            elif page == 4:
                model[row][0] = '---'

        set_changed(gui, page, True)
        if page == 2:
            gui.iface.view_qal.grab_focus()
            gtk_qal.start_trim(gui)
        elif page == 4:
            gui.iface.view_gbl.grab_focus()
            gtk_gbl.drop_seqs(gui)


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
