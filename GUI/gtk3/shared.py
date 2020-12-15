# 2020 Leo Kaindl

import hashlib
import logging
from argparse import Namespace
from pathlib import Path

import gi

from ab12phylo import msa
from gtk_msa import PAGE

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject

from GUI.gtk3 import gtk_io, gtk_rgx, gtk_qal, gtk_msa, gtk_gbl, gtk_ml

LOG = logging.getLogger(__name__)
TOOLS = Path(__file__).resolve().parents[2] / 'ab12phylo' / 'tools'
USER = 'leo.kaindl@tum.de'
SEP = '?*?'
RAW_MSA = Path('Trim') / 'raw_msa.fasta'
IMPORT_MSA = Path('import') / 'import_raw_msa.fasta'
MSA = 'msa.fasta'
MISSING = 'missing_samples.tsv'
PREVIEW = Path('Trim') / 'trim_preview.png'
CBAR = Path('Trim') / 'colorbar.png'
LEFT = Path('Trim') / 'msa_gbl_pre.png'
RIGHT = Path('Trim') / 'msa_gbl_post.png'
ALPHA = .25
DPI = 300
BUF_SIZE = 128 * 1024
NUCLEOTIDES = ['A', 'C', 'G', 'T', 'N', 'else', '-', ' ', 'S']
KXLIN = {
    'A': (0.92, 1, 0.4, 1),
    'C': (0.46, 1, 0.44, 1),
    'G': (0.16, 0.44, 0.8, 1),
    'T': (1, 0.47, 0.66, 1),
    'N': (0.84, 0.84, 0.84, 0.6),
    'else': (0, 0, 0, 1),
    '-': (1, 1, 1, 0),
    ' ': (1, 1, 1, 0),
    'S': (1, 1, 1, 0)}

# TODO continue project variables
# re-fresh page content. automatically called
REFRESH = [gtk_io.refresh, gtk_rgx.refresh, gtk_qal.refresh, gtk_msa.refresh, gtk_gbl.refresh, gtk_ml.refresh]
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
        # gui.iface.next.get_style_context().remove_class(Gtk.STYLE_CLASS_SUGGESTED_ACTION)
    else:
        # enable the next page
        iface.notebook.get_children()[page + 1].set_sensitive(True)
        # set earlier pages to 'unchanged'
        data.change_indicator[:page + 1] = [False] * (page + 1)
        # gui.iface.next.get_style_context().remove_class(Gtk.STYLE_CLASS_SUGGESTED_ACTION)


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


def get_hashes(gui, msa_path):
    msa_hash = file_hash(gui.wd / msa_path)
    if msa_hash != gui.data.msa_hash:
        set_changed(gui, PAGE)
        gui.data.msa_hash = msa_hash


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
    # REFRESH[data.page](gui)
    LOG.debug('stepped back to page %d' % data.page)


def re_run(gui, *args):
    """
    Depending on the currently visible page, re-run the matching background task.
    Handles the Refresh button.
    """
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    # gui.iface.refresh.grab_focus()
    RERUN[page](gui)


def refresh(gui, *args):
    """
    Call the refresh function of the current page and hide or show
    the appropriate buttons.
    :param gui:
    :param args:
    :return:
    """
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    REFRESH[page](gui)
    # hide or show these two actions depending on applicability
    gui.iface.refresh.props.visible = bool(page in RERUN)
    gui.iface.gene_roll.props.visible = bool(page in SELECT)


def select_gene_and_redo(gui, *args):
    """
    A page-independent handler for selecting a different gene.
    Currently, iface.gene_roll is only visible from one page, so a bit useless.
    :param gui:
    :param args:
    :return:
    """
    if args:
        gui = args[-1]
    page = gui.iface.notebook.get_current_page()
    if page == 2:  # trim preview
        gtk_qal.parse(gui.iface.gene_roll, None, gui)
        gtk_qal.start_trim(gui)
    # TODO continue for new pages


def get_column(list_store, col_idx):
    """Extract a column from a Gtk.ListStore, because they are annoying."""
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col


def load_image(zoom_ns, page, gtk_bin, img_path, width, height):
    """
    Load an image into a GtkImage child inside gtk_bin, keeping track of observed sizes.
    :param zoom_ns: where all the zooming info is stored, specifically [min_w, min_h,
    w_now, h_now] of the image inside :param gtk_bin inside the zoom_ns.sizes dict.
    """
    # save the current size: only called on enlarging -> easy
    zoom_ns.sizes[page] = [width, height, width, height]

    child = gtk_bin.get_child()
    if type(child) != Gtk.Image:
        gtk_bin.remove(child)
        child = Gtk.Image()
        gtk_bin.add(child)

    pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(img_path), width=max(1, width), height=height,
        preserve_aspect_ratio=False)
    child.set_from_pixbuf(pb)


def load_colorbar(gtk_image, wd):
    gtk_image.set_from_pixbuf(GdkPixbuf.Pixbuf.new_from_file_at_scale(
        str(wd / CBAR), width=250, height=100, preserve_aspect_ratio=True))


def get_height_resize(widget, event, spacer, scroll_wins, lower=0):
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
    return h


def get_hadj(iface):
    return iface.zoom.adj.get_value() * 2


def scale(gtk_bin, x, y):
    """
    Stretch a pixbuf residing in a GtkImage, the single child of
    :param gtk_bin:
    :param x: horizontal scaling factor
    :param y: vertical scaling factor
    :return:
    """
    child = gtk_bin.get_child()
    if type(child) != Gtk.Image:
        LOG.info('won\'t rescale %s' % str(type(child)))
        return False
    pb = child.get_pixbuf()
    child.set_from_pixbuf(pb.scale_simple(
        min(14000, pb.get_width() * x), min(14000, pb.get_height() * y),
        GdkPixbuf.InterpType.BILINEAR))
    return True


def x_scale(adj, gui, zoom_ns):
    """Horizontally scale a preview."""
    with GObject.signal_handler_block(adj, zoom_ns.handle):
        a = adj.props
        page = gui.iface.notebook.get_current_page()
        data, iface = gui.data, gui.iface

        a.value = max(.2, a.value)
        x = a.value / zoom_ns.bak
        min_w, min_h, w_now, h_now = zoom_ns.sizes[page]

        if x * w_now > 2 * min_w:
            LOG.debug('re-loading images')
            # load larger, so zoom in won't happen so soon again
            a.value = min(a.upper, a.value + 2 * a.step_increment)
            x = a.value / zoom_ns.bak
            if page == 2:
                load_image(zoom_ns, page, iface.qal_eventbox, gui.wd / PREVIEW,
                           data.qal_shape[0] * a.value * 2, h_now)
            elif page == 4:
                load_image(zoom_ns, page, iface.gbl_left_vp, gui.wd / LEFT,
                           data.msa_shape[0] * a.value * 2, h_now)
                load_image(zoom_ns, page, iface.gbl_right_vp, gui.wd / RIGHT,
                           data.msa_shape[2] * a.value * 2, h_now)
            zoom_ns.sizes[page] = [w_now * x, min_h, w_now * x, h_now]
        else:
            LOG.debug('scale x: %.2f fold' % x)
            if page == 2:
                scale(iface.qal_eventbox, x, 1)
            elif page == 4:
                scale(iface.gbl_left_vp, x, 1)
                scale(iface.gbl_right_vp, x, 1)
            zoom_ns.sizes[page] = [min(w_now * x, min_w), min_h, w_now * x, h_now]
        zoom_ns.bak = a.value


def xy_scale(widget, event, gui, page):
    """
    Handles zoom in / zoom out on Ctrl+mouse wheel in both x and y
    :param widget:
    :param event:
    :param gui:
    :param page:
    :return:
    """
    data, iface = gui.data, gui.iface
    accel_mask = Gtk.accelerator_get_default_mod_mask()
    if event.state & accel_mask == Gdk.ModifierType.CONTROL_MASK:
        with GObject.signal_handler_block(iface.zoom.adj, iface.zoom.handle):
            direction = event.get_scroll_deltas()[2]
            a = iface.zoom.adj.props
            bak = a.value
            min_w, min_h, w_now, h_now = iface.zoom.sizes[page]
            if direction > 0:  # scrolling down -> zoom out -> simple
                # go one tick down
                a.value = max(.2, a.value - a.step_increment)
                if a.value == bak:
                    return
                new = a.value / bak
                LOG.debug('scale xy: %.2f fold, %.1f' % (new, a.value))
                if page == 2:
                    scale(iface.qal_eventbox, new, new)
                elif page == 4:
                    scale(iface.gbl_left_vp, new, new)
                    scale(iface.gbl_right_vp, new, new)
                # adjust the saved sizes
                iface.zoom.sizes[page] = [min(min_w * new, min_w), min(min_h * new, min_h), w_now * new, h_now * new]
            else:
                a.value = min(a.upper, a.value + a.step_increment)
                new = a.value / bak
                if a.value == bak:
                    return
                min_w, min_h, w_now, h_now = iface.zoom.sizes[page]
                if w_now / min_w * new > 3 or h_now / min_h * new > 4:
                    # re-load is due. jump ahead by incrementing again
                    a.value = min(a.upper, a.value + 2 * a.step_increment)
                    new = a.value / bak
                    LOG.debug('re-loading images, scale xy: %.2f fold, %.1f' % (new, a.value))
                    if page == 2:
                        load_image(iface.zoom, page, iface.qal_eventbox, gui.wd / PREVIEW,
                                   data.qal_shape[0] * a.value * 2, data.qal_shape[1])
                    elif page == 4:
                        load_image(iface.zoom, page, iface.gbl_left_vp, gui.wd / LEFT,
                                   data.msa_shape[0] * a.value * 2, data.gbl_shape[1])
                        load_image(iface.zoom, page, iface.gbl_right_vp, gui.wd / RIGHT,
                                   data.msa_shape[2] * a.value * 2, data.gbl_shape[1])
                else:
                    LOG.debug('scale xy: %.2f fold, %.1f' % (new, a.value))
                    # scale the easy way
                    if page == 2:
                        scale(iface.qal_eventbox, new, new)
                    elif page == 4:
                        scale(iface.gbl_left_vp, new, new)
                        scale(iface.gbl_right_vp, new, new)
                iface.zoom.sizes[page] = [min_w, min_h, w_now * new, h_now * new]
            iface.zoom.bak = a.value


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


def select_seqs(event_box, loc, page, zoom_ns, tv, ns):
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
    if tv.props.visible:
        tv.grab_focus()
    sel = tv.get_selection()

    accel_mask = Gtk.accelerator_get_default_mod_mask()
    rect, baseline = event_box.get_allocated_size()
    mo, tree_path_iterator = sel.get_selected_rows()
    idcs = {tp[0] for tp in tree_path_iterator}

    h_now = zoom_ns.sizes[page][3]
    idx = int((loc.y - (rect.height - h_now) / 2) / h_now * len(mo))
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


def init_gene_roll(gui):
    """
    Initialize gene switcher combo box with the previously selected gene from the dataset.
    """
    data, iface = gui.data, gui.iface
    with iface.gene_roll.handler_block(iface.gene_handler):
        iface.gene_roll.remove_all()
        genes = list(data.genes)
        [iface.gene_roll.append_text(gene) for gene in genes]
        if len(genes) > 1:
            iface.gene_roll.insert_text(0, 'all')
            genes.insert(0, 'all')
        if data.gene_for_preview:
            idx = genes.index(data.gene_for_preview)
        else:
            idx = 0
        iface.gene_roll.set_active(idx)


def keep_visible(sel, adj, ns):
    """
    For keyboard navigation in previews, scroll the TreeView and keep the selection up-to-date
    :param sel: a TreeSelection
    :param adj: the Props of a Vadjustment in a ScrolledWindow
    :param ns: a Namespace where the previous selection state is stored (amongst others)
    :return:
    """
    mo, tp_iter = sel.get_selected_rows()
    tps = {tp[0] for tp in tp_iter}
    if 'sel' not in ns:
        ns.sel = tps
        return

    tp = tps - ns.sel
    ns.sel = tps
    # scrolling only for incremental keyboard selection
    if len(tp) != 1:
        return
    tp = tp.pop()

    if (tp + 1) / len(mo) > (adj.value + adj.page_size) / adj.upper:
        # scroll down
        adj.value = min(adj.upper - adj.page_size,
                        ((tp + 1) / (len(mo) + 1) - adj.page_size / adj.upper) * adj.upper)
    elif tp / len(mo) < adj.value / adj.upper:
        # scroll up
        adj.value = tp / (len(mo) + 1) * adj.upper


def get_cmd(algo, gui, remote=False):
    """
    Initializes a commandline msa_build object
    :param algo:
    :param gui:
    :param remote:
    :return:
    """
    data, iface = gui.data, gui.iface
    args = Namespace(**{
        'dir': gui.wd,
        'genes': data.genes,
        'msa_algo': algo,
        'user': USER,
        'msa': gui.wd / MSA,
        'sep': SEP,
        'missing_samples': None
    })
    iface.aligner = msa.msa_build(args, None, no_run=True)
    if remote:
        cmd = iface.aligner.build_remote('%s', no_run=True)
    else:
        cmd = iface.aligner.build_local('%s', no_run=True)
    return cmd
