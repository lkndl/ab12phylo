import gi
import logging
from pathlib import Path

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk, Gdk as gdk

from GUI.gtk3 import commons

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 0
FILETYPES = ['.ab1', '.seq', '.fasta']


def init(data, iface):
    # trace file types
    # create a TreeView model
    data.file_type_model = gtk.ListStore(str, bool)
    [data.file_type_model.append([filetype, False]) for filetype in FILETYPES]

    # check ABI traces by default
    data.file_type_model[0][1] = True
    data.filetypes.add('.ab1')

    iface.view_filetypes.set_model(data.file_type_model)
    iface.view_filetypes.set_headers_visible(False)
    iface.view_filetypes.append_column(
        gtk.TreeViewColumn(title='Filetype', cell_renderer=gtk.CellRendererText(), text=0))
    crt = gtk.CellRendererToggle()
    crt.connect('toggled', change_filetypes, data)
    iface.view_filetypes.append_column(
        gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))

    for mo, tv, file_type, file_types in zip([data.trace_model, data.csv_model],
                                             [iface.view_trace_path, iface.view_csv_path],
                                             ['trace', 'csv'], [data.filetypes, {'.csv'}]):
        # data.__setattr__(model, gtk.ListStore(str))
        # mo, tv = data.__getattribute__(model), iface.__getattribute__(tree_view)
        tv.set_model(mo)
        tv.set_headers_visible(False)
        tv.append_column(gtk.TreeViewColumn(title='Paths',
                                            cell_renderer=gtk.CellRendererText(), text=0))

        sel = tv.get_selection()
        sel.set_mode(gtk.SelectionMode.MULTIPLE)

        iface.__getattribute__('add_%s_folder' % file_type) \
            .connect('clicked', add_folder, iface, data, file_types, mo)
        iface.__getattribute__('add_%s_manual' % file_type) \
            .connect('clicked', add_manually, iface, data, mo)
        try:
            iface.__getattribute__('add_%s_whitelist' % file_type) \
                .connect('clicked', add_manually, iface, data, mo, True)
        except AttributeError:
            pass  # wellsplate whitelist not planned

        iface.__getattribute__('remove_%s' % file_type) \
            .connect('clicked', commons.delete_rows, iface, data, PAGE, sel)
        iface.__getattribute__('delete_all_%s' % file_type) \
            .connect('clicked', commons.delete_rows, iface, data, PAGE, sel, mo)

    iface.files_next.connect('clicked', commons.proceed, data, iface)


def add_folder(widget, iface, data, file_types, model):
    dialog = gtk.FileChooserDialog(title='select folder(s)', parent=None, select_multiple=True,
                                   action=gtk.FileChooserAction.SELECT_FOLDER)
    dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                       gtk.STOCK_OPEN, gtk.ResponseType.OK)
    duplicates, new_paths = list(), list()
    response = dialog.run()
    if response == gtk.ResponseType.OK:
        for folder in dialog.get_filenames():
            try:
                folder = Path(folder).resolve()
                # get all the matching files
                new_paths = [list(folder.glob('**/*%s' % ext)) for ext in file_types]
                # flatten the LoL
                new_paths = [str(item) for sublist in new_paths for item in sublist]
                # append to the ListStore
                [model.append([path]) for path in new_paths if [path] not in model]
            except UnicodeDecodeError as ex:
                LOG.info(ex)
        commons.refresh_files(iface, data, PAGE)
        LOG.debug('found %d new paths in folder(s)' % len(new_paths))

    elif response == gtk.ResponseType.CANCEL:
        LOG.debug('cancel')
    elif response == gtk.ResponseType.DELETE_EVENT:
        LOG.debug('escape')

    dialog.destroy()

    if duplicates:
        commons.show_message_dialog('Some files already selected', duplicates)


def add_manually(widget, iface, data, model, whitelists=False):
    dialog = gtk.FileChooserDialog(title='select files', parent=None,
                                   action=gtk.FileChooserAction.OPEN, select_multiple=True)
    dialog.add_buttons(gtk.STOCK_CANCEL, gtk.ResponseType.CANCEL,
                       gtk.STOCK_OPEN, gtk.ResponseType.OK)
    not_found, duplicates = list(), list()
    response = dialog.run()
    if response == gtk.ResponseType.OK:
        # read the filenames
        new_paths = dialog.get_filenames()
        if whitelists:
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
                    new_paths.append(str(path))
                else:
                    not_found.append(string_path)

        # append to the ListStore
        [model.append([path]) for path in new_paths if [path] not in model]
        commons.refresh_files(iface, data, PAGE)
        LOG.info('added %d paths' % len(new_paths))

    elif response == gtk.ResponseType.CANCEL:
        LOG.info('cancel')
    elif response == gtk.ResponseType.DELETE_EVENT:
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


def add_new_entries(store, entries_to_add, iface):
    """
    Adds new entries to a list and returns the modified list
    as well as lists of both newly added entries and observed duplicates
    :param store: the trace_paths list of a dataset object
    :param entries_to_add: obvs
    :param iface: the namespace containing all named widgets of a gui object
    :return:
    """
    dups, news = list(), list()
    for entry in entries_to_add:
        if entry in store:
            dups.append(entry)
        else:
            store.append(entry)
            news.append(entry)
    if len(news) > 0:
        commons.set_changed(iface, PAGE, True)
    return store, news, dups
