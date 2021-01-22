# 2020 Leo Kaindl

import hashlib
import logging
import re
import threading
from argparse import Namespace
from time import sleep

import gi
from time import time
import pandas as pd

from ab12phylo import msa
from ab12phylo_gui.static import PATHS, USER, SEP, BUF_SIZE

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject

from ab12phylo_gui import gtk_io, gtk_rgx, gtk_qal, gtk_msa, gtk_gbl, gtk_blast, gtk_ml, gtk_tree

LOG = logging.getLogger(__name__)

# re-fresh page content. automatically called
REFRESH = [module.refresh for module in [gtk_io, gtk_rgx, gtk_qal, gtk_msa,
                                         gtk_gbl, gtk_blast, gtk_ml, gtk_tree]]
# re-run background threads. -> "REFRESH" button
RERUN = {1: gtk_rgx.start_read, 2: gtk_qal.start_trim,
         4: gtk_gbl.start_gbl, 7: gtk_tree.start_phy}
# where the gene selector is visible
SELECT = {2}

regex = re.compile('\\.[\\d]+$')


def _inc_priv_timestamp():
    return str(time()).replace('.', '')[6:]


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
    nb = iface.notebook
    if changed:
        # disable later pages
        [page.set_sensitive(False) for page in nb.get_children()[page + 1: nb.get_n_pages()]]
        # set later pages to 'changed'
        data.change_indicator[page:] = [True] * (nb.get_n_pages() - page)
        data.change_indicator[:page] = [False] * page
    else:
        # enable the next page
        nb.get_children()[page + 1].set_sensitive(True)
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


def show_notification(gui, msg, items=None, stay_secs=False):
    revealer, iface = gui.iface.revealer, gui.iface
    iface.reveal_title.set_text(msg)
    iface.reveal_list.props.parent.props.visible = items
    if items:
        iface.reveal_list.props.buffer.props.text = '\n'.join(items)
    gui.iface.revealer.set_transition_duration(stay_secs / 4 if stay_secs else 250)
    revealer.set_reveal_child(True)
    if stay_secs:
        threading.Thread(target=hide, args=[revealer, stay_secs]).start()


def hide(revealer, stay_for):
    sleep(stay_for)
    revealer.set_reveal_child(False)


def tv_keypress(widget, event, gui, page, selection):
    if Gdk.keyval_name(event.keyval) == 'Delete':
        LOG.debug('delete selected row')
        delete_files_from_input_selection(widget, gui, page, selection)
    else:
        LOG.debug('registered keypress, did nothig')


def get_hashes(gui, msa_path, page):
    msa_hash = file_hash(gui.wd / msa_path)
    if msa_hash != gui.data.msa_hash:
        set_changed(gui, page)
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
    :param widget: optional, for use as callback
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
            gtk_io.refresh(gui)
            gtk_rgx.reset_columns(gui, do_parse=True)
        elif page == 1:
            if 1 < sum(data.rx_fired) < 5:
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
        elif page == 5:
            iface.tempspace.df.to_csv(
                gui.wd / PATHS.tsv, sep='\t', na_rep='', header=True, index=True)
            set_changed(gui, 5, False)
        elif page == 6:
            gtk_tree.start_phy(gui)
        set_changed(gui, page, False)

    # then proceed
    iface.notebook.next_page()
    data.page = iface.notebook.get_current_page()
    LOG.debug('proceeded to page %d' % data.page)


def step_back(widget, gui):
    """Go back to the previous page; set it to un-changed; re-view"""
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


def get_column(list_store, col_idx):
    """Extract a column from a Gtk.ListStore, because they are annoying"""
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col


def load_image(zoom_ns, page, gtk_bin, img_path, w=None, h=None):
    """
    Load an image into a GtkImage child inside gtk_bin, keeping track of observed sizes.
    :param zoom_ns: where all the zooming info is stored, specifically [min_w, min_h,
    w_now, h_now] of the image inside :param gtk_bin inside the zoom_ns.sizes dict.
    """
    # save the current size: only called on enlarging -> easy
    zoom_ns.sizes[page] = [w, h, w, h]

    child = gtk_bin.get_child()
    if child:
        gtk_bin.remove(child)
    child = Gtk.Image()
    gtk_bin.add(child)
    path = str(img_path)

    if h and w:
        pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            filename=path, width=max(1, w), height=h,
            preserve_aspect_ratio=False)
    elif w:
        pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            path, width=w, height=-1, preserve_aspect_ratio=True)
    elif h:
        pb = GdkPixbuf.Pixbuf.new_from_file_at_scale(
            filename=str(path), width=-1, height=h, preserve_aspect_ratio=True)
    else:
        raise ValueError('Specify at least one out of width / height')
    child.set_from_pixbuf(pb)


def load_colorbar(gtk_image, wd, gbar=False):
    try:
        path = PATHS.gbar if gbar else PATHS.cbar
        gtk_image.set_from_pixbuf(GdkPixbuf.Pixbuf.new_from_file_at_scale(
            str(wd / path), width=250, height=100, preserve_aspect_ratio=True))
    except FileNotFoundError:
        pass


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
    return iface.zoomer.adj.get_value() * 2


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
    """Horizontally scale a preview"""
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
                load_image(zoom_ns, page, iface.qal_eventbox, gui.wd / PATHS.preview,
                           data.qal_shape[0] * a.value * 2, h_now)
            elif page == 4:
                load_image(zoom_ns, page, iface.gbl_left_vp, gui.wd / PATHS.left,
                           data.msa_shape[0] * a.value * 2, h_now)
                load_image(zoom_ns, page, iface.gbl_right_vp, gui.wd / PATHS.right,
                           data.msa_shape[2] * a.value * 2, h_now)
            elif page == 7:
                load_image(zoom_ns, page, iface.msa_eventbox, gui.wd / PATHS.phylo_msa,
                           data.phy.shape[0] * a.value * 2, h_now)
            zoom_ns.sizes[page] = [w_now * x, min_h, w_now * x, h_now]
        else:
            LOG.debug('scale x: %.2f fold' % x)
            if page == 2:
                scale(iface.qal_eventbox, x, 1)
            elif page == 4:
                scale(iface.gbl_left_vp, x, 1)
                scale(iface.gbl_right_vp, x, 1)
            elif page == 7:
                scale(iface.msa_eventbox, x, 1)
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
        with GObject.signal_handler_block(iface.zoomer.adj, iface.zoomer.handle):
            direction = event.get_scroll_deltas()[2]
            a = iface.zoomer.adj.props
            bak = a.value
            min_w, min_h, w_now, h_now = iface.zoomer.sizes[page]
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
                elif page == 7:
                    scale(iface.msa_eventbox, new, new)
                # adjust the saved sizes
                iface.zoomer.sizes[page] = [min(min_w * new, min_w), min(min_h * new, min_h), w_now * new, h_now * new]
            else:
                a.value = min(a.upper, a.value + a.step_increment)
                new = a.value / bak
                if a.value == bak:
                    return
                min_w, min_h, w_now, h_now = iface.zoomer.sizes[page]
                if w_now / min_w * new > 3 or h_now / min_h * new > 4:
                    # re-load is due. jump ahead by incrementing again
                    a.value = min(a.upper, a.value + 2 * a.step_increment)
                    new = a.value / bak
                    LOG.debug('re-loading images, scale xy: %.2f fold, %.1f' % (new, a.value))
                    if page == 2:
                        load_image(iface.zoomer, page, iface.qal_eventbox, gui.wd / PATHS.preview,
                                   data.qal_shape[0] * a.value * 2, data.qal_shape[1])
                    elif page == 4:
                        load_image(iface.zoomer, page, iface.gbl_left_vp, gui.wd / PATHS.left,
                                   data.msa_shape[0] * a.value * 2, data.gbl_shape[1])
                        load_image(iface.zoomer, page, iface.gbl_right_vp, gui.wd / PATHS.right,
                                   data.msa_shape[2] * a.value * 2, data.gbl_shape[1])
                    elif page == 7:
                        load_image(iface.zoomer, page, iface.msa_eventbox, gui.wd / PATHS.phylo_msa,
                                   data.phy.shape[0] * a.value * 2, data.phy.shape[1])
                else:
                    LOG.debug('scale xy: %.2f fold, %.1f' % (new, a.value))
                    # scale the easy way
                    if page == 2:
                        scale(iface.qal_eventbox, new, new)
                    elif page == 4:
                        scale(iface.gbl_left_vp, new, new)
                        scale(iface.gbl_right_vp, new, new)
                    elif page == 7:
                        scale(iface.msa_eventbox, new, new)
                iface.zoomer.sizes[page] = [min_w, min_h, w_now * new, h_now * new]
            if page == 7:
                iface.tree_pane.set_position(300)
            iface.zoomer.bak = a.value


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
    if iface.thread.is_alive():
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


def update_ML(iface, page, ml):
    """
    Extract the numbers of finished tree searches / bootstrap
    iterations from the RAxML output to update the progress bar.
    """
    if ml.key:
        seen_set = ml.seen[ml.key]
        motif = ml.motifs[ml.key]
        for line in ml.stdout[len(ml.stdout) - 20:]:
            if motif in line:
                try:
                    seen_set.add(int(line.split(motif)[1].split(',')[0]))
                except ValueError:
                    pass
        iface.i = len(seen_set) + ml.prev
    return update(iface, page)


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
    # LOG.debug('select_seqs')
    if tv.props.visible:
        tv.grab_focus()
    sel = tv.get_selection()

    accel_mask = Gtk.accelerator_get_default_mod_mask()
    rect, baseline = event_box.get_allocated_size()
    mo, tree_path_iterator = sel.get_selected_rows()
    idcs = {tp[0] for tp in tree_path_iterator}

    h_now = zoom_ns.sizes[page][3]
    # print('%d:%d:%d:%d' % (tv.get_margin_bottom(), rect.height, h_now, loc.y))
    idx = int((loc.y - (rect.height - h_now) / 2) / (h_now - tv.get_margin_bottom()) * len(mo))
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


def save_row_edits(cell, path, new_text, tv, col):
    LOG.debug('save_row_edits')
    mo = tv.get_model()
    old_text = mo[path][col]
    if old_text == new_text:
        return
    mo[path][col] = new_text


def delete_files_from_input_selection(widget, gui, page, selection, delete_all=False):
    if delete_all:
        delete_all.clear()
        if widget == gui.iface.delete_all_trace:
            gui.data.trace_store.clear()
    else:
        model, iterator = selection.get_selected_rows()
        [model.remove(model.get_iter(row)) for row in reversed(sorted(iterator))]
    set_changed(gui, page, True)
    REFRESH[page](gui)


def delete_and_ignore_rows(widget, event, gui, page, sel, ns):
    """
    Keep track of the rows that will not be written to the next fasta and delete them from the treeview.
    """
    if Gdk.keyval_name(event.keyval) == 'Delete':
        LOG.debug('delete_and_ignore_rows')
        model, tree_path_iterator = sel.get_selected_rows()
        if 'ignore_ids' not in ns:
            ns.ignore_ids = set() if page == 4 else dict()
        for row in reversed(sorted(tree_path_iterator)):
            if page == 2:
                _id, gene = model[row][:2]
                _id = shift_versions_down_on_deletion(
                    _id, gui.data.seqdata[gene], gui.data.metadata[gene])
                if gene in ns.ignore_ids:
                    ns.ignore_ids[gene].add(_id)
                else:
                    ns.ignore_ids[gene] = {_id}
                model.remove(model.get_iter(row))
            elif page == 4:
                ns.ignore_ids.add(model[row][0])
                model[row][0] = '---'

        set_changed(gui, page, True)
        if page == 2:
            gui.iface.view_qal.grab_focus()
            gtk_qal.start_trim(gui)
        elif page == 4:
            gui.iface.view_gbl.grab_focus()
            gtk_gbl.drop_seqs(gui)


def shift_versions_down_on_deletion(_id, genedata, genemeta):
    match = regex.search(_id)
    if match:
        next_v = int(_id[match.start() + 1:]) + 1
        stem = _id[:match.start()]
    else:
        next_v = 1
        stem = _id

    # move entry out of the way
    r = genedata.pop(_id)
    r.id = _id + '_' + _inc_priv_timestamp()
    genedata[r.id] = r
    genemeta[r.id] = genemeta.pop(_id)
    genemeta[r.id]['quality'] = 'manually dropped at trim_data'

    next_id = '%s.%d' % (stem, next_v)
    while next_id in genedata:
        # re-organize seqdata and metadata
        r = genedata.pop(next_id)
        r.id = stem if next_v == 1 else '%s.%d' % (stem, next_v - 1)
        genedata[r.id] = r
        genemeta[r.id] = genemeta.pop(next_id)

        next_v += 1
        next_id = '%s.%d' % (stem, next_v)

    return '%s.%d' % (stem, next_v - 1) if next_v > 1 else _id


def write_metadata(gui):
    """
    Write the metadata dictionary to the metadata.tsv
    """
    md = gui.data.metadata
    # convert metadata dict do pandas DataFrame
    df = pd.concat({gene: pd.DataFrame.from_dict(
        md[gene], orient='index') for gene in md.keys()})
    df.index.names = ['gene', 'id']
    df.reset_index(level=0, inplace=True)
    # write to file
    df.to_csv(gui.wd / PATHS.tsv, sep='\t', na_rep='', header=True, index=True)


def init_gene_roll(gui):
    """
    Initialize gene switcher combo box with the previously selected gene from the dataset.
    """
    data, iface = gui.data, gui.iface
    with iface.gene_roll.handler_block(iface.gene_handler):
        iface.gene_roll.remove_all()
        genes = list(data.genes)
        if not genes:
            return
        [iface.gene_roll.append_text(gene) for gene in genes]
        if len(genes) > 1:
            iface.gene_roll.insert_text(0, 'all')
            genes.insert(0, 'all')
        if data.gene_for_preview and data.gene_for_preview in genes:
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


def get_msa_build_cmd(algo, wd, genes, remote=False):
    """
    Initializes a commandline msa_build object
    :param algo:
    :param gui:
    :param remote:
    :return:
    """
    args = Namespace(**{
        'dir': wd,
        'genes': genes,
        'msa_algo': algo,
        'user': USER,
        'msa': wd / PATHS.msa,
        'sep': SEP,
        'missing_samples': None
    })
    aligner = msa.msa_build(args, None, no_run=True)
    if remote:
        cmd = aligner.build_remote('%s', no_run=True)
    else:
        cmd = aligner.build_local('%s', no_run=True)
    return aligner, cmd


class bump_log_level:

    def __init__(self, log, off=True):
        self.level = max(log.level, logging.INFO)
        self.off = off

    def __enter__(self):
        if not self.off:
            logging.disable(self.level)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        logging.disable(logging.NOTSET)


def edit_numerical_entry_up_down(widget, event):
    key = Gdk.keyval_name(event.keyval)
    if key == 'Up':
        widget.set_text(str(1 + int(widget.get_text())))
        return True
    elif key == 'Down':
        widget.set_text(str(max(0, -1 + int(widget.get_text()))))
        return True


def edit_numerical_entry(widget):
    # filter for numbers only
    value = ''.join([c for c in widget.get_text() if c.isdigit()])
    widget.set_text(value if value else '')  # entering 0 is annoying
