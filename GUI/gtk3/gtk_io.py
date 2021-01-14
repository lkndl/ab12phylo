# 2020 Leo Kaindl

import logging
from argparse import Namespace
from pathlib import Path

import gi

from GUI.gtk3 import gtk_proj, shared

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 0
FILETYPES = ['.ab1', '.seq', '.fasta', '.fa']


def init(gui):
    data, iface = gui.data, gui.iface

    # blink the progress bar so the window allocates a size
    iface.prog_bar.props.visible = True
    gui.win.show_all()
    iface.prog_bar.props.visible = False

    # trace file types
    # create a TreeView model
    iface.file_type_model = gtk_proj.picklable_liststore(str, bool)
    [iface.file_type_model.append([file_type, False]) for file_type in FILETYPES]

    # check ABI traces by default
    iface.file_type_model[0][1] = True

    iface.view_filetypes.set_model(iface.file_type_model)
    iface.view_filetypes.set_headers_visible(False)
    iface.view_filetypes.append_column(
        Gtk.TreeViewColumn(title='Filetype', cell_renderer=Gtk.CellRendererText(), text=0))
    crt = Gtk.CellRendererToggle()
    crt.connect('toggled', lambda widget, path: iface.file_type_model.set(
        iface.file_type_model.get_iter(path), [1], [not iface.file_type_model[path][1]]))
    iface.view_filetypes.append_column(
        Gtk.TreeViewColumn(title='Selected', cell_renderer=crt, active=1))

    # trace files
    iface.add_refs.connect('clicked', add_manually, gui, data.trace_store, 'references')

    # initial or old number of entries in model used for scrolling to end after adding
    iface.file_nums = dict()

    # set up the file paths tables
    for mo, tv, file_type, col in zip([data.trace_store, data.plate_store],
                                      [iface.view_trace_path, iface.view_csv_path],
                                      ['trace', 'csv'], [7, 3]):
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
            .connect('clicked', add_folder, gui, file_type, mo)
        iface.__getattribute__('add_%s_manual' % file_type) \
            .connect('clicked', add_manually, gui, mo)
        try:
            iface.__getattribute__('add_%s_whitelist' % file_type) \
                .connect('clicked', add_manually, gui, mo, 'whitelists')
        except AttributeError:
            pass  # wellsplate whitelist not planned

        iface.__getattribute__('remove_%s' % file_type) \
            .connect('clicked', shared.delete_files_from_input_selection, gui, PAGE, sel)
        tv.connect('key-press-event', shared.tv_keypress, gui, PAGE, sel)
        iface.__getattribute__('delete_all_%s' % file_type) \
            .connect('clicked', shared.delete_files_from_input_selection, gui, PAGE, sel, mo)

    iface.rasterize.set_visible(False)
    # horizontal scaling
    iface.zoom = Namespace()
    iface.zoom.adj = iface.qal_scale.get_adjustment()
    iface.zoom.adj.configure(1, .2, 2.6, .2, 0, 0)
    iface.qal_scale.set_digits(1)
    iface.gbl_scale.set_digits(1)
    iface.qal_scale.add_mark(1, Gtk.PositionType.BOTTOM, None)
    iface.gbl_scale.add_mark(1, Gtk.PositionType.BOTTOM, None)
    iface.zoom.bak = iface.zoom.adj.get_value()
    iface.zoom.handle = iface.zoom.adj.connect_after('value-changed', shared.x_scale, gui, iface.zoom)
    iface.zoom.sizes = dict()  # will save with the image parent name here


def add_folder(widget, gui, file_type, model):
    data, iface = gui.data, gui.iface

    if file_type == 'trace':
        file_types = {a[0] for a in shared.get_column(iface.file_type_model, (0, 1)) if a[1]}
    else:
        file_types = {'.csv'}

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
                model, new_paths, duplicates = add_new_entries(model, new_paths, gui, kw)
            except UnicodeDecodeError as ex:
                LOG.info(ex)
        refresh(gui)
        LOG.debug('found %d new paths in folder(s)' % len(new_paths))
    dialog.destroy()

    if duplicates:
        shared.show_notification(gui, 'Files already selected:', duplicates)


def add_manually(widget, gui, model, *args):
    """
    Load trace, wellsplate or reference paths into the appropriate GtkListStore,
    also remembering if this is a reference.
    """
    data, iface = gui.data, gui.iface
    args = ['nope'] if not args else args
    dialog = Gtk.FileChooserDialog(title='select files', parent=None,
                                   action=Gtk.FileChooserAction.OPEN, select_multiple=True)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                       Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
    not_found, duplicates = list(), list()
    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        # read the filenames
        filenames = dialog.get_filenames()
        if 'whitelists' in args:
            # read lines from the files
            extracted = list()
            for whitelist in filenames:
                try:
                    extracted += open(whitelist, 'r').read().strip().split('\n')
                except UnicodeDecodeError as ex:
                    LOG.info(ex)
            filenames = extracted

        # sort extracted into valid and invalid file paths
        new_paths = list()
        for string_path in filenames:
            path = Path(string_path).resolve()
            if path.is_file():
                new_paths.append(path)
            else:
                not_found.append(string_path)

        # where to write
        if widget in [iface.add_csv_folder, iface.add_csv_manual]:
            kw = 'csv'
        elif widget == iface.add_refs:
            kw = 'ref'
        else:
            kw = 'trace'
        # append to ListStore
        model, new_paths, duplicates = add_new_entries(model, new_paths, gui, kw)
        refresh(gui)
        LOG.info('added %d paths' % len(new_paths))
    dialog.destroy()

    if not_found or duplicates:
        shared.show_notification(gui, 'File troubles', ['not found:%s' % f for f in not_found] +
                                 ['duplicate:%s' % f for f in duplicates])


def add_new_entries(model, new_paths, gui, *args):
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
    data, iface = gui.data, gui.iface
    if 'ref' in args:
        color = iface.AQUA
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
        shared.set_changed(gui, PAGE)
    return model, news, dups


def scroll_to_end(widget, rectangle, iface, tv, mo):
    """After new entries have been added to it, the TreeView will scroll to its end."""
    if iface.notebook.get_current_page() != 0:
        return
    if mo not in iface.file_nums or len(mo) > iface.file_nums[mo]:
        adj = tv.get_vadjustment()
        adj.set_value(adj.get_upper() - adj.get_page_size())
    iface.file_nums[mo] = len(mo)


def refresh(gui, page=PAGE):
    """
    Adjusts the display of the number of selected trace files on the files page
    and toggles reading wellsplates or not, inactivating respective GUI elements.
    :return:
    """
    data, iface = gui.data, gui.iface
    num_traces = len(data.trace_store)
    if num_traces > 0:
        iface.trace_number.set_label('%d files' % num_traces)
        iface.trace_number.set_visible(True)
    else:
        iface.trace_number.set_label('0 files')
        iface.trace_number.set_visible(False)
