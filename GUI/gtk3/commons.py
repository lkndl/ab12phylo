import gi
from pathlib import Path
import logging

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from GUI.gtk3 import regex, quality

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)


def get_changed(iface, page):
    return iface.change_indicator[page]


def get_errors(iface, page):
    return iface.errors_indicator[page]


def set_changed(iface, page, changed):
    """
    Set page to changed and disable further stages; set page to unchanged and enable next page
    :param iface: the namespace containing all named widgets of a gui object
    :param page: the index of the current page
    :param changed: toggle
    :return:
    """
    if changed:
        # disable later pages
        [page.set_sensitive(False) for page in iface.notebook.get_children()[page + 1: iface.notebook.get_n_pages()]]
        # set later pages to 'changed'
        iface.change_indicator[page:] = [True] * (iface.notebook.get_n_pages() - page)
    else:
        # enable the next page
        iface.notebook.get_children()[page + 1].set_sensitive(True)
        # set earlier pages to 'unchanged'
        iface.change_indicator[:page + 1] = [False] * (page + 1)


def set_errors(iface, page, errors):
    iface.errors_indicator[page] = errors


def show_message_dialog(message, list_to_print=None):
    dialog = Gtk.MessageDialog(transient_for=None, flags=0,
                               buttons=Gtk.ButtonsType.OK, message_type=Gtk.MessageType.WARNING,
                               text=message)
    # dialog.format_secondary_text('\n'.join(not_found))  # not copyable -> not user-friendly
    if list_to_print:
        txt_buf = Gtk.TextBuffer()
        txt_buf.set_text('\n'.join(list_to_print))
        dialog.get_message_area().add(Gtk.TextView().new_with_buffer(txt_buf))
    dialog.show_all()  # important
    dialog.run()
    dialog.destroy()
    LOG.debug('showed a message')


def refresh_files(iface, data, page):
    """
    Adjusts the display of the number of selected trace files on the files page
    and toggles reading wellsplates or not, inactivating respective GUI elements.
    :return:
    """
    num_traces = len(data.trace_model)
    if num_traces > 0:
        iface.trace_number.set_label('%d files' % num_traces)
        iface.trace_number.set_visible(True)
    else:
        iface.trace_number.set_label('0 files')
        iface.trace_number.set_visible(False)

    # toggle _reading_plates
    iface.plates = len(data.csv_model) > 0
    # self._set_page_changed(not iface.plates or self._get_page_changed())
    set_changed(iface, page, not iface.plates or get_changed(iface, page))
    [iface.__getattribute__(name).set_sensitive(iface.plates)
     for name in ['plate_regex_label', 'plate_rx', 'wp_rx_desc', 'wp_lbl',
                  'wp_rx', 'wellsplate_buttons', 'wellsplate_regex_box']]
    # toggle radiobutton line back off
    if not iface.triple_rt.get_active():
        [iface.__getattribute__(name).set_sensitive(False) for name in
         ['plate_rx', 'plate_regex_label']]


def delete_rows(widget, iface, data, page, selection, delete_all=False):
    if delete_all:
        delete_all.clear()
    else:
        model, iterator = selection.get_selected_rows()
        [model.remove(model.get_iter(row)) for row in reversed(sorted(iterator))]
    set_changed(iface, page, True)
    refresh_files(iface, data, page)


def proceed(widget, gui):
    data, iface = gui.data, gui.interface
    page = iface.notebook.get_current_page()
    # first integrate changes to the dataset
    if get_changed(iface, page):
        if page == 0:
            regex.reset(gui)
        # TODO
        elif page == 1:
            # check if everything seems fine on the
            if get_errors(iface, page):
                show_message_dialog('There are still errors on the page!')
                return
            quality.reset(gui)
        set_changed(iface, page, False)

    # then proceed
    iface.notebook.next_page()
    LOG.debug('proceeded to page %d' % iface.notebook.get_current_page())


def step_back(widget, gui):
    data, iface = gui.data, gui.interface
    iface.notebook.prev_page()
    page = iface.notebook.get_current_page()
    set_changed(iface, page, False)
    LOG.debug('stepped back to page %d' % page)


def get_column(list_store, col_idx):
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col
