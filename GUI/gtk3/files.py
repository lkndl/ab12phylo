import logging
from pathlib import Path

import gi

from GUI.gtk3 import commons

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 0
FILETYPES = ['.ab1', '.seq', '.fasta', '.fa']


def init(gui):
    data, iface = gui.data, gui.interface

    # MARK trace file types
    # create a TreeView model
    data.file_type_model = Gtk.ListStore(str, bool)
    [data.file_type_model.append([file_type, False]) for file_type in FILETYPES]

    # check ABI traces by default
    data.file_type_model[0][1] = True
    data.filetypes.add('.ab1')

    iface.view_filetypes.set_model(data.file_type_model)
    iface.view_filetypes.set_headers_visible(False)
    iface.view_filetypes.append_column(
        Gtk.TreeViewColumn(title='Filetype', cell_renderer=Gtk.CellRendererText(), text=0))
    crt = Gtk.CellRendererToggle()
    crt.connect('toggled', change_filetypes, data)
    iface.view_filetypes.append_column(
        Gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))

    # MARK trace files
    iface.add_refs.connect('clicked', add_manually, gui, data.trace_store, 'references')

    # initial or old number of entries in model used for scrolling to end after adding
    iface.file_nums = dict()

    # set up the file paths tables
    for mo, tv, file_type, file_types, col in zip([data.trace_store, data.plate_store],
                                                  [iface.view_trace_path, iface.view_csv_path],
                                                  ['trace', 'csv'], [data.filetypes, {'.csv'}], [7, 3]):
        tv.set_model(mo)
        tv.set_headers_visible(False)
        tv.append_column(Gtk.TreeViewColumn(title='Paths',
                                            cell_renderer=Gtk.CellRendererText(),
                                            foreground=col, text=0))
        tv.connect('size-allocate', scroll_to_end, iface, tv, mo)
        iface.file_nums[mo] = 0
        sel = tv.get_selection()
        sel.set_mode(Gtk.SelectionMode.MULTIPLE)

        iface.__getattribute__('add_%s_folder' % file_type) \
            .connect('clicked', add_folder, gui, file_types, mo)
        iface.__getattribute__('add_%s_manual' % file_type) \
            .connect('clicked', add_manually, gui, mo)
        try:
            iface.__getattribute__('add_%s_whitelist' % file_type) \
                .connect('clicked', add_manually, gui, mo, 'whitelists')
        except AttributeError:
            pass  # wellsplate whitelist not planned

        iface.__getattribute__('remove_%s' % file_type) \
            .connect('clicked', commons.delete_rows, gui, PAGE, sel)
        iface.__getattribute__('delete_all_%s' % file_type) \
            .connect('clicked', commons.delete_rows, gui, PAGE, sel, mo)

    iface.files_next.connect('clicked', commons.proceed, gui)


def add_folder(widget, gui, file_types, model):
    data, iface = gui.data, gui.interface
    dialog = Gtk.FileChooserDialog(title='select folder(s)', parent=None, select_multiple=True,
                                   action=Gtk.FileChooserAction.SELECT_FOLDER)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                       Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
    duplicates, new_paths = list(), list()
    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        for folder in dialog.get_filenames():
            try:
                folder = Path(folder).resolve()
                # get all the matching files
                new_paths = [list(folder.glob('**/*%s' % ext)) for ext in file_types]
                # flatten the LoL
                new_paths = [item for sublist in new_paths for item in sublist]
                # where to write
                kw = 'csv' if '.csv' in file_types else 'trace'
                # append to the ListStore
                model, new_paths, duplicates = add_new_entries(model, new_paths, iface, kw)
            except UnicodeDecodeError as ex:
                LOG.info(ex)
        commons.refresh_files(gui, PAGE)
        LOG.debug('found %d new paths in folder(s)' % len(new_paths))

    elif response == Gtk.ResponseType.CANCEL:
        LOG.debug('cancel')
    elif response == Gtk.ResponseType.DELETE_EVENT:
        LOG.debug('escape')

    dialog.destroy()

    if duplicates:
        commons.show_message_dialog('Some files already selected', duplicates)


def add_manually(widget, gui, model, *args):
    """
    Load trace or wellsplate paths into the appropriate GtkListStore,
    also remembering if this is a reference.
    """
    data, iface = gui.data, gui.interface
    args = ['nope'] if not args else args
    dialog = Gtk.FileChooserDialog(title='select files', parent=None,
                                   action=Gtk.FileChooserAction.OPEN, select_multiple=True)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                       Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
    not_found, duplicates = list(), list()
    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        # read the filenames
        new_paths = dialog.get_filenames()
        if 'whitelists' in args:
            # read lines from the files
            extracted = list()
            for whitelist in new_paths:
                try:
                    extracted += open(whitelist, 'r').read().strip().split('\n')
                except UnicodeDecodeError as ex:
                    LOG.info(ex)

            # sort extracted into valid and invalid file paths

            new_paths = list()
            for string_path in extracted:
                path = Path(string_path).resolve()
                if path.is_file():
                    new_paths.append(path)
                else:
                    not_found.append(string_path)
        else:
            new_paths = [Path(string_path) for string_path in new_paths]

        # where to write
        if widget in [iface.add_csv_folder, iface.add_csv_manual]:
            kw = 'csv'
        elif widget == iface.add_refs:
            kw = 'ref'
        else:
            kw = 'trace'
        # append to ListStore
        model, new_paths, duplicates = add_new_entries(model, new_paths, iface, kw)
        commons.refresh_files(gui, PAGE)
        LOG.info('added %d paths' % len(new_paths))

    elif response == Gtk.ResponseType.CANCEL:
        LOG.info('cancel')
    elif response == Gtk.ResponseType.DELETE_EVENT:
        LOG.info('escape')

    dialog.destroy()

    if not_found:
        commons.show_message_dialog('Some file paths were invalid', not_found)
    if duplicates:
        commons.show_message_dialog('Some files already selected', duplicates)


def change_filetypes(widget, path, data):
    # toggle the button
    data.file_type_model[path][1] = not data.file_type_model[path][1]
    # adapt the filetypes to read for the next folder
    if data.file_type_model[path][1]:
        data.filetypes.add(data.file_type_model[path][0])
    else:
        data.filetypes.remove(data.file_type_model[path][0])
    LOG.debug('selected ' + ' and '.join(data.filetypes))


def add_new_entries(model, new_paths, iface, *args):
    """
    Adds new entries to a list and returns the modified list
    as well as lists of both newly added entries and observed duplicates.
    If *args, the new entries are reference sequences and the value in
    the color column will change accordingly
    :param model: the ListStore
    :param new_paths: obvs
    :param iface: the namespace containing all named widgets of a gui object
    :return:
    """
    if 'ref' in args:
        color = iface.BLUE
        is_ref = True
    elif 'trace' in args or 'csv' in args:
        color = iface.FG
        is_ref = False
    else:
        assert False

    old_entries = [line[0] for line in model]
    dups, news = list(), list()

    for ppath in new_paths:
        path = str(ppath)
        if path in old_entries:
            dups.append(path)
        else:
            if 'csv' in args:
                model.append([path, ppath.name, '', color])
            else:  # trace file
                model.append([path, ppath.name, '', '', '', is_ref, False, color])
            news.append(path)
    if len(news) > 0:
        commons.set_changed(iface, PAGE, True)
    return model, news, dups


def scroll_to_end(widget, rectangle, iface, tv, mo):
    """After new entries have been added to it, the TreeView will scroll to its end."""
    if len(mo) > iface.file_nums[mo]:
        LOG.debug('scrolling to end)')
        adj = tv.get_vadjustment()
        adj.set_value(adj.get_upper() - adj.get_page_size())
    iface.file_nums[mo] = len(mo)
