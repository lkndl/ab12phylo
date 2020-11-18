# 2020 Leo Kaindl

import gi
from pathlib import Path
import logging

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from GUI.gtk3 import files, regex, quality

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)


def get_changed(gui, page):
    return gui.data.change_indicator[page]


def get_errors(gui, page):
    return gui.data.errors_indicator[page]


def set_errors(gui, page, errors):
    gui.data.errors_indicator[page] = errors


def set_changed(gui, page, changed=True):
    """
    Set page to changed and disable further stages; set page to unchanged and enable next page
    :param gui:
    :param page: the index of the current page
    :param changed: toggle
    :return:
    """
    data, iface = gui.data, gui.interface
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
        txt_buf.set_text('\n'.join(list_to_print))
        dialog.get_message_area().add(Gtk.TextView().new_with_buffer(txt_buf))
    dialog.add_button('OK', 0)
    dialog.get_widget_for_response(0).grab_focus()
    dialog.show_all()  # important
    dialog.run()
    dialog.destroy()
    LOG.debug('showed a message')


def delete_rows(widget, gui, page, selection, delete_all=False):
    data, iface = gui.data, gui.interface
    if delete_all:
        delete_all.clear()
    else:
        model, iterator = selection.get_selected_rows()
        [model.remove(model.get_iter(row)) for row in reversed(sorted(iterator))]
    set_changed(gui, page, True)
    files.refresh_files(gui, page)


def proceed(widget, gui):
    data, iface = gui.data, gui.interface
    page = iface.notebook.get_current_page()
    # first integrate changes to the dataset
    if get_changed(gui, page):
        if page == 0:
            regex.reset(gui, do_parse=True)
        elif page == 1:
            # check if everything ok
            regex.re_check(gui)
            if get_errors(gui, page):
                show_message_dialog('There are still errors on the page!')
                return
            if 1 < sum(iface.rx_fired) < 5:
                print(iface.rx_fired)
                show_message_dialog('Make sure all columns have been parsed.')
                return
            regex.read_files(gui)
            return  # leave this alone
        elif page == 2:
            quality.trim_all(gui)
            # TODO
            print('wohooo')
        set_changed(gui, page, False)

    # then proceed
    iface.notebook.next_page()
    LOG.debug('proceeded to page %d' % iface.notebook.get_current_page())


def step_back(widget, gui):
    data, iface = gui.data, gui.interface
    iface.notebook.prev_page()
    page = iface.notebook.get_current_page()
    set_changed(gui, page, False)
    LOG.debug('stepped back to page %d' % page)


def get_column(list_store, col_idx):
    col = list()
    for row in list_store:
        col.append(row[col_idx])
    return col


def update(iface, bar, page):
    """
    Keeps a progress bar up-to-date
    :param iface:
    :param bar: the progress bar to update
    :param page: the page to freeze
    :return:
    """

    if iface.running:
        iface.notebook.get_children()[page].set_sensitive(False)
        bar.set_visible(True)
        bar.set_fraction(iface.frac)
        bar.set_text(iface.txt)
        return True
    else:
        iface.notebook.get_children()[page].set_sensitive(True)
        bar.set_visible(False)
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
        except Exception as ex:
            LOG.exception(ex)


def _unpickle_liststore(cls, col_types, rows):
    inst = cls.__new__(cls)
    inst.__init__(*col_types)
    for row in rows:
        inst.append(row)
    return inst
