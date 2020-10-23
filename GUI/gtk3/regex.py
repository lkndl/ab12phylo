import json
import logging
import webbrowser
from pathlib import Path
from time import sleep

import gi
import re
import requests

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from GUI.gtk3 import commons

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 1
ERRORS = ['<span foreground="blue">groups?</span>',
          '<span foreground="red">no match</span>',
          '<span foreground="green">use groups!</span>']
MARKUP = ['<span foreground="blue">reverse</span>',
          '<span foreground="red">no match</span>',
          '',
          '<span foreground="green">use groups!</span>']


def init(gui):
    data, iface = gui.data, gui.interface
    iface.single_rt.connect('toggled', rx_toggle, data, iface)
    iface.triple_rt.connect('toggled', rx_toggle, data, iface)
    iface.triple_rt.join_group(iface.single_rt)
    iface.single_rt.set_active(True)

    iface.regex_apply.connect('clicked', parse_triple, data, iface)
    iface.wp_apply.connect('clicked', parse_single, data, iface, iface.wp_rx, 0)
    iface.single_rx.connect('activate', parse_triple, data, iface)

    for thing, col in zip(['well', 'gene', 'plate', 'wp'], [0, 2, 1, 0]):
        entry = iface.__getattribute__('%s_rx' % thing)
        entry.connect('activate', parse_single, data, iface, entry, col)

    iface.reverse_rx_chk.connect('toggled', rev_adjust, data, iface)
    iface.reverse_rx.connect('activate', parse_single, data, iface, iface.reverse_rx, 3)
    iface.regex_try_online.connect('clicked', try_online, data, iface)

    # allow deletion
    for widget, sel in zip([iface.remove_path_regex, iface.remove_csv_regex],
                           [iface.view_trace_regex.get_selection,
                            iface.view_csv_regex.get_selection]):
        sel().set_mode(Gtk.SelectionMode.MULTIPLE)
        widget.connect('clicked', commons.delete_rows, iface, data, PAGE, sel())

    # connect buttons
    iface.regex_next.connect('clicked', commons.proceed, gui)
    iface.regex_back.connect('clicked', commons.step_back, gui)
    reset(gui)
    commons.refresh_files(iface, data, PAGE)


def reset(gui):
    data, iface = gui.data, gui.interface
    # create TreeView model
    # path extraction: well/id, (plate), gene, (reversed), file, path
    data.rx_model = Gtk.ListStore(str, str, str, str, str, str)
    # fill last column with initial filename data
    [data.rx_model.append([''] * 4 + [Path(path[0]).name] + [path[0]]) for path in data.trace_model]
    if iface.plates:
        # plate ID, filename, path
        data.wp_model = Gtk.ListStore(str, str, str)
        [data.wp_model.append([''] + [Path(path[0]).name] + [path[0]]) for path in data.csv_model]
    # iface.view_trace_regex.set_headers_visible(False)
    iface.view_trace_regex.set_model(data.rx_model)
    iface.view_csv_regex.set_model(data.wp_model)
    # remove old columns:
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        [tree_view.remove_column(col) for col in tree_view.get_columns()]

    # columns depend on _reading_plates
    if iface.plates:
        for title, column in zip(['well', 'plate', 'gene', 'file'], [0, 1, 2, 4]):
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(), markup=column))
        # wellsplates:
        for title, column in zip(['plate ID', 'file'], list(range(2))):
            iface.view_csv_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(), markup=column))
        iface.rx_fired = False, False
    else:
        for title, column in zip(['sample', 'gene', 'file'], [0, 2, 4]):
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(), markup=column))
        iface.rx_fired = False, True

    # reset_sort_size(data, iface)
    iface.reverse_rx_chk.set_active(False)
    iface.search_rev = False


def parse_single(widget, data, iface, entry, col):
    LOG.debug('parsing from %s' % entry.get_name())
    regex = re.compile(entry.get_text())
    errors = False or commons.get_errors(iface, PAGE)
    changed = False or commons.get_changed(iface, PAGE)

    if widget in [iface.wp_rx, iface.wp_apply]:
        model = data.wp_model
        iface.rx_fired = True, iface.rx_fired[1]
    else:
        model = data.rx_model

    if entry is not iface.reverse_rx:
        for i, row in enumerate(model):
            file = row[-2]
            try:
                m = regex.search(file).groups()[0]
                if not changed:
                    changed = m != row[col]
                model[i][col] = m
            except ValueError as ve:
                # maybe wrong number of groups
                model[i][col] = ERRORS[0]
                errors, changed = True, True
            except AttributeError as ae:
                # no match
                model[i][col] = ERRORS[1]
                errors, changed = True, True
            except IndexError as ie:
                # no groups used
                model[i][col] = ERRORS[2]
                errors, changed = True, True
    else:
        for i, row in enumerate(model):
            file = row[-2]
            try:
                m = regex.search(file).groups()[0]
                # reverse read
                if not changed:
                    changed = MARKUP[0] != row[col]
                model[i][col] = MARKUP[0]
            except ValueError as ve:
                # no match
                errors, changed = True, True
                model[i][col] = MARKUP[1]
            except AttributeError as ae:
                # forward read
                if not changed:
                    changed = MARKUP[2] != row[col]
                model[i][col] = MARKUP[2]
            except IndexError as ie:
                # wrong groups
                errors, changed = True, True
                model[i][col] = MARKUP[3]
    commons.set_changed(iface, PAGE, changed)
    commons.set_errors(iface, PAGE, errors)


def parse_triple(widget, data, iface):
    if iface.single_rt.get_active():
        LOG.debug('parsing with single regex')
        # if parsing with only a single regex
        regex = re.compile(iface.single_rx.get_text())
        errors = False or commons.get_errors(iface, PAGE)
        changed = False or commons.get_changed(iface, PAGE)

        for idx, row in enumerate(data.rx_model):
            file = row[-2]
            try:
                m = regex.search(file)
                plate, gene, well = m.groups() if iface.plates else (None, *m.groups())
                if not changed:
                    changed = not bool(row == [well, plate, gene, row[3], file])
                data.rx_model[idx] = [well, plate, gene, row[3], file, row[-1]]
            except ValueError as ve:
                errors, changed = True, True
                data.rx_model[idx][0] = ERRORS[0]
            except AttributeError as ae:
                errors, changed = True, True
                data.rx_model[idx][0] = ERRORS[1]
            except IndexError as ie:
                errors, changed = True, True
                data.rx_model[idx][0] = ERRORS[2]
        commons.set_changed(iface, PAGE, changed)
        commons.set_errors(iface, PAGE, errors)
    else:
        parse_single(None, data, iface, iface.well_rx, 0)
        parse_single(None, data, iface, iface.gene_rx, 2)
        if iface.plates:
            parse_single(None, data, iface, iface.plate_rx, 1)
    if iface.search_rev:
        parse_single(None, data, iface, iface.reverse_rx, 3)
    iface.rx_fired = True, iface.rx_fired[1]
    reset_sort_size(data, iface)
    LOG.debug('parse_triple done')


def rx_toggle(widget, data, iface):
    # (In)activates the Entry fields and causes a parse event
    if widget.get_active():
        three = ['plate_rx', 'gene_rx', 'well_rx', 'plate_regex_label', 'gene_regex_label', 'well_regex_label']
        if widget == iface.single_rt:
            iface.single_rx.set_sensitive(True)
            [iface.__getattribute__(widget).set_sensitive(False) for widget in three]
        else:
            iface.single_rx.set_sensitive(False)
            [iface.__getattribute__(widget).set_sensitive(True) for widget in three]
        # adjust the sensitivity of the plate regex line
        if not iface.plates:
            iface.plate_rx.set_sensitive(False)
            iface.plate_regex_label.set_sensitive(False)
        LOG.debug('toggled radiobutton %s' % widget.get_name())
        # parse again
        parse_triple(widget, data, iface)


def rev_adjust(widget, data, iface):
    n_cols = iface.view_trace_regex.get_n_columns()
    if widget.get_active():
        iface.search_rev = True
        # enable entry field
        iface.reverse_rx.set_sensitive(True)
        # create new column
        # px = Gtk.CellRendererPixbuf()
        # col = Gtk.TreeViewColumn(title=' ', cell_renderer=px)
        # col.set_cell_data_func(px,
        #                        lambda column, cell, model, _iter, user_data:
        #                        cell.set_property('icon-name', model[_iter][3]))
        # iface.view_trace_regex.insert_column(col, iface.view_trace_regex.get_n_columns() - 1)
        iface.view_trace_regex.insert_column(
            Gtk.TreeViewColumn(title='reverse', cell_renderer=Gtk.CellRendererText(), markup=3), n_cols - 1)
        # cause parsing
        parse_single(None, data, iface, iface.reverse_rx, 3)
    else:
        iface.search_rev = False
        iface.reverse_rx.set_sensitive(False)
        for col in iface.view_trace_regex.get_columns():
            if col.get_title() == 'reverse':
                iface.view_trace_regex.remove_column(col)

    reset_sort_size(data, iface)


def reset_sort_size(data, iface):
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        tree_view.columns_autosize()
        for col_index in range(tree_view.get_n_columns()):
            tree_view.get_column(col_index).set_sort_column_id(col_index)

    # check the dataset for empty strings
    if '' not in commons.get_column(data.rx_model, 0) \
            + commons.get_column(data.rx_model, 2) \
            + commons.get_column(data.wp_model, 0):
        if not iface.plates or '' not in commons.get_column(data.rx_model, 1):
            commons.set_errors(iface, PAGE, False)
            LOG.debug('found no errors')
            return
    assert commons.get_errors(iface, PAGE) or not iface.rx_fired[0] or not iface.rx_fired[1]


def try_online(widget, data, iface):
    LOG.debug('generating online help')
    url = 'https://regex101.com/api/regex/'
    headers = {'content-type': 'application/json'}

    payload = dict(flavor='python', flags='gm', delimiter='"')
    payload['regex'] = iface.single_rx.get_text()
    content = tuple(i.get_text() for i in [iface.wp_rx, iface.single_rx, iface.well_rx,
                                           iface.gene_rx, iface.plate_rx, iface.reverse_rx]) \
              + ('\n'.join([row[-2] for row in data.wp_model]),
                 '\n'.join([row[-2] for row in data.rx_model]))

    payload['testString'] = 'wellsplate ID:\n%s\n\nsingle RegEx:\n%s\n\n' \
                            'separate expressions:\n%s\n%s\n%s\n\n' \
                            'reverse read:\n%s\n\n' \
                            'plate filenames:\n%s\n\n' \
                            'trace filenames:\n%s' % content

    payload = json.dumps(payload, indent=2)
    try:
        response = requests.request("POST", url, data=payload, headers=headers)
        result = json.loads(response.text)

        webbrowser.open('https://regex101.com/r/' + result['permalinkFragment'])

        # wait for the page to load
        sleep(4)
        # delete the entry
        payload = json.dumps(dict(deleteCode=result['deleteCode']), indent=2)
        response = requests.request("DELETE", url, data=payload, headers=headers)
    except Exception:
        commons.show_message_dialog('online help failed')
