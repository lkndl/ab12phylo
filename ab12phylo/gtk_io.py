# 2021 Leo Kaindl

import logging
from pathlib import Path

import gi
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from ab12phylo.gtk_base import ab12phylo_app_base
from ab12phylo.gtk_proj import picklable_liststore

LOG = logging.getLogger(__name__)
PAGE = 0


class io_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        # blink the progress bar so the window allocates a size
        iface.prog_bar.props.visible = True
        iface.win.show_all()
        iface.prog_bar.props.visible = False

        # trace file types
        # create a TreeView model
        iface.file_type_model = picklable_liststore(str, bool)
        [iface.file_type_model.append([file_type, False]) for
         file_type in ['.ab1', '.seq', '.fasta', '.fa', '.txt']]

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
        iface.add_refs.connect('clicked', self.add_manually, data.trace_store, 'references')

        # initial or old number of entries in model used for scrolling to end after adding
        iface.file_nums = dict()

        # set up the file paths tables
        for mo, tv, file_type, col in zip([data.trace_store, data.plate_store],
                                          [iface.view_trace_path, iface.view_csv_path],
                                          ['trace', 'csv'], [7, 3]):
            tv.set_model(mo)
            tv.set_headers_visible(False)
            tv.append_column(Gtk.TreeViewColumn(title='Paths',
                                                cell_renderer=Gtk.CellRendererText(editable=True),
                                                foreground=col, text=0))
            tv.connect('size-allocate', self.scroll_to_end, tv, mo)
            iface.file_nums[mo] = 0
            sel = tv.get_selection()
            sel.set_mode(Gtk.SelectionMode.MULTIPLE)

            iface.__getattribute__('add_%s_folder' % file_type) \
                .connect('clicked', self.add_folder, file_type, mo)
            iface.__getattribute__('add_%s_manual' % file_type) \
                .connect('clicked', self.add_manually, mo)
            try:
                iface.__getattribute__('add_%s_whitelist' % file_type) \
                    .connect('clicked', self.add_manually, mo, 'whitelists')
            except AttributeError:
                pass  # wellsplate whitelist not planned

            iface.__getattribute__('remove_%s' % file_type) \
                .connect('clicked', self.delete_files_from_input_selection, PAGE, sel)
            tv.connect('key-press-event', self.tv_keypress, PAGE, sel)
            iface.__getattribute__('delete_all_%s' % file_type) \
                .connect('clicked', self.delete_files_from_input_selection, PAGE, sel, mo)

        # horizontal scaling
        iface.zoomer.adj = iface.qal_scale.get_adjustment()
        iface.zoomer.adj.configure(1, .2, 2.6, .2, 0, 0)
        iface.qal_scale.set_digits(1)
        iface.gbl_scale.set_digits(1)
        iface.tree_scale.set_digits(1)
        iface.qal_scale.add_mark(1, Gtk.PositionType.BOTTOM, None)
        iface.gbl_scale.add_mark(1, Gtk.PositionType.BOTTOM, None)
        iface.tree_scale.add_mark(1, Gtk.PositionType.BOTTOM, None)
        iface.zoomer.bak = iface.zoomer.adj.get_value()
        iface.zoomer.handle = iface.zoomer.adj.connect_after('value-changed', self.x_scale, iface.zoomer)
        iface.zoomer.sizes = dict()  # will save with the image parent name here

    def add_folder(self, widget, file_type, model):
        if file_type == 'trace':
            file_types = {a[0] for a in self.iface.
                file_type_model.get_column((0, 1)) if a[1]}
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
                    model, new_paths, duplicates = self.add_new_entries(model, new_paths, kw)
                except UnicodeDecodeError as ex:
                    LOG.info(ex)
            self.refresh()
            LOG.debug('found %d new paths in folder(s)' % len(new_paths))
        dialog.destroy()

        if duplicates:
            self.show_notification('Files already selected:', duplicates)

    def add_manually(self, widget, model, *args):
        """
        Load trace, wellsplate or reference paths into the appropriate GtkListStore,
        also remembering if this is a reference.
        """
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
            if widget in [self.iface.add_csv_folder, self.iface.add_csv_manual]:
                kw = 'csv'
            elif widget == self.iface.add_refs:
                kw = 'ref'
            else:
                kw = 'trace'
            # append to ListStore
            model, new_paths, duplicates = self.add_new_entries(model, new_paths, kw)
            self.refresh()
            LOG.info('added %d paths' % len(new_paths))
        dialog.destroy()

        if not_found or duplicates:
            self.show_notification('File troubles', ['not found:%s' % f for f in not_found] +
                                   ['duplicate:%s' % f for f in duplicates], secs=0)

    def add_new_entries(self, model, new_paths, *args):
        """
        Adds new entries to a list and returns the modified list
        as well as lists of both newly added entries and observed duplicates.
        If *args, the new entries are reference sequences and the value in
        the color column will change accordingly
        :param model: the ListStore
        :param new_paths: obvs
        :param iface: the namespace containing all named widgets of a ab12phylo_app_base object
        :return:
        """
        if 'ref' in args:
            color = self.iface.AQUA
            is_ref = True
        elif 'trace' in args or 'csv' in args:
            color = None  # iface.FG
            is_ref = False
        else:
            raise RuntimeWarning(f'weird file type, neither ref, trace or csv')

        old_entries = [line[0] for line in model]
        dups, news = list(), list()

        for ppath in new_paths:
            path = str(ppath)
            if path in old_entries:
                dups.append(path)
            else:
                if 'csv' in args:
                    model.append([path, ppath.name, '', color])
                elif ppath.suffix.lower() != '.ab1' and not is_ref:
                    # no ABI trace, no ref and several records in the file -> list separately
                    rids = [r.id for r in SeqIO.parse(ppath, 'fasta')]
                    if len(rids) > 1:
                        [model.append([path + '~~' + rid, rid, '', '', '',
                                       False, False, color]) for rid in rids]
                    else:
                        model.append([path, ppath.name, '', '', '', is_ref, False, color])
                else:  # ABI trace file
                    model.append([path, ppath.name, '', '', '', is_ref, False, color])
                news.append(path)
        if len(news) > 0:
            self.set_changed(PAGE)
        return model, news, dups

    def scroll_to_end(self, widget, rectangle, tv, mo):
        """After new entries have been added to it, the TreeView will scroll to its end"""
        if self.iface.notebook.get_current_page() != 0:
            return
        if mo not in self.iface.file_nums or len(mo) > self.iface.file_nums[mo]:
            adj = tv.get_vadjustment()
            adj.set_value(adj.get_upper() - adj.get_page_size())
        self.iface.file_nums[mo] = len(mo)

    def refresh(self):
        """
        Adjusts the display of the number of selected trace files on the files page
        and toggles reading wellsplates or not, inactivating respective GUI elements.
        :return:
        """
        num_traces = len(self.data.trace_store)
        if num_traces > 0:
            self.iface.trace_number.set_label('%d files' % num_traces)
            self.iface.trace_number.set_visible(True)
        else:
            self.iface.trace_number.set_label('0 files')
            self.iface.trace_number.set_visible(False)
