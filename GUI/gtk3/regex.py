import json
import logging
import webbrowser
from pathlib import Path
from time import sleep

import gi
import re
import requests

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

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


# TODO errors on page after removing lines
# TODO try for gene in references


def init(gui):
    data, iface = gui.data, gui.interface
    iface.plates, iface.search_rev = True, False

    iface.single_rt.connect('toggled', rx_toggle, gui)
    iface.triple_rt.connect('toggled', rx_toggle, gui)
    iface.triple_rt.join_group(iface.single_rt)
    iface.single_rt.set_active(True)

    iface.regex_apply.connect('clicked', parse_triple, gui)
    iface.wp_apply.connect('clicked', parse_single, gui, iface.wp_rx, 2)
    iface.single_rx.connect('activate', parse_triple, gui)

    for thing, col in zip(['well', 'gene', 'plate', 'wp'], [2, 4, 3, 2]):
        entry = iface.__getattribute__('%s_rx' % thing)
        entry.connect('activate', parse_single, gui, entry, col)
        entry.connect('focus_out_event', parse_single, gui, entry, col)

    iface.reverse_rx_chk.connect('toggled', rev_adjust, gui)
    iface.reverse_rx.connect('activate', parse_single, gui, iface.reverse_rx, 6)
    iface.regex_try_online.connect('clicked', try_online, gui)

    # allow deletion
    for widget, sel in zip([iface.remove_path_regex, iface.remove_csv_regex],
                           [iface.view_trace_regex.get_selection,
                            iface.view_csv_regex.get_selection]):
        sel().set_mode(Gtk.SelectionMode.MULTIPLE)
        widget.connect('clicked', commons.delete_rows, gui, PAGE, sel())

    # connect buttons
    iface.regex_next.connect('clicked', commons.proceed, gui)
    iface.regex_back.connect('clicked', commons.step_back, gui)
    reset(gui)
    commons.refresh_files(gui, PAGE)


def reset(gui):
    data, iface = gui.data, gui.interface
    iface.view_trace_regex.set_model(data.trace_store)
    iface.view_csv_regex.set_model(data.plate_store)
    # remove old columns:
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        [tree_view.remove_column(col) for col in tree_view.get_columns()]

    # columns depend on _reading_plates
    if iface.plates:
        for title, column in zip(['well', 'plate', 'gene', 'file'], [2, 3, 4, 1]):
            crm = Gtk.CellRendererText(editable=True)
            crm.connect('edited', cell_edit, iface.view_trace_regex, column)
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crm, text=column,
                                   foreground_rgba=7))
        # wellsplates:
        for title, column in zip(['plate ID', 'file'], [2, 1]):
            crm = Gtk.CellRendererText(editable=True)
            crm.connect('edited', cell_edit, iface.view_csv_regex, column)
            iface.view_csv_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crm, text=column))
        iface.rx_fired = False, False
    else:
        for title, column in zip(['sample', 'gene', 'file'], [2, 4, 1]):
            crm = Gtk.CellRendererText(editable=True)
            crm.connect('edited', cell_edit, iface.view_trace_regex, column)
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crm, text=column,
                                   foreground_rgba=7))
        iface.rx_fired = False, True

    # reset_sort_size(gui)
    iface.reverse_rx_chk.set_active(False)
    iface.search_rev = False

    # fire the initial parse
    parse_single(iface.wp_rx, gui, iface.wp_rx, 2)
    parse_triple(None, gui)


def parse_single(widget, gui, entry, col):
    # TODO reverse
    data, iface = gui.data, gui.interface
    LOG.debug('parsing from %s' % entry.get_name())
    regex = re.compile(entry.get_text())
    errors = False or commons.get_errors(iface, PAGE)
    changed = False or commons.get_changed(iface, PAGE)

    if widget in [iface.wp_rx, iface.wp_apply]:
        model = data.plate_store
        iface.rx_fired = True, iface.rx_fired[1]
        traces = False
    else:
        model = data.trace_store
        traces = True

    if entry is not iface.reverse_rx:
        for i, row in enumerate(model):
            # skip references
            if traces and row[5]:
                continue
            file = row[1]
            try:
                m = regex.search(file).groups()[0]
                if not changed:
                    changed = m != row[col]
                model[i][col] = m
                continue
            except ValueError as ve:
                # maybe wrong number of groups
                model[i][col] = 'wrong number of groups?'
            except AttributeError as ae:
                # no match
                model[i][col] = 'no match'
            except IndexError as ie:
                # no groups used
                'use groups'
            model[i][-1] = iface.RED
            errors, changed = True, True
    else:
        for i, row in enumerate(model):
            # skip references
            if traces and row[5]:
                continue
            file = row[1]
            try:
                m = regex.search(file).groups()[0]
                # reverse read
                if not changed:
                    changed = row[col] is not True
                model[i][col] = True
            except AttributeError as ae:
                # forward read
                if not changed:
                    changed = row[col] is not False
                model[i][col] = False
            except (IndexError, ValueError):
                model[i][-1] = iface.RED
                errors, changed = True, True
                model[i][col] = ''
            except Exception:
                assert False
    commons.set_changed(iface, PAGE, changed)
    commons.set_errors(iface, PAGE, errors)


def parse_triple(widget, gui):
    data, iface = gui.data, gui.interface
    if iface.single_rt.get_active():
        LOG.debug('parsing with single regex')
        # if parsing with only a single regex
        regex = re.compile(iface.single_rx.get_text())
        errors = False or commons.get_errors(iface, PAGE)
        changed = False or commons.get_changed(iface, PAGE)

        for idx, row in enumerate(data.trace_store):
            if row[4]:
                continue
            file = row[-2]
            try:
                m = regex.search(file)
                plate, gene, well = m.groups() if iface.plates else (None, *m.groups())
                if not changed:
                    # TODO continue here
                    changed = not bool(row == [well, plate, gene, row[3], file])
                data.trace_store[idx] = row[:2] + [well, plate, gene] + row[5:]
                # data.trace_store[idx] = [well, plate, gene] + row[3:]
            except ValueError as ve:
                errors, changed = True, True
                data.trace_store[idx][0] = ERRORS[0]
            except AttributeError as ae:
                errors, changed = True, True
                data.trace_store[idx][0] = ERRORS[1]
            except IndexError as ie:
                errors, changed = True, True
                data.trace_store[idx][0] = ERRORS[2]
        commons.set_changed(iface, PAGE, changed)
        commons.set_errors(iface, PAGE, errors)
    else:
        parse_single(None, gui, iface.well_rx, 0)
        parse_single(None, gui, iface.gene_rx, 2)
        if iface.plates:
            parse_single(None, gui, iface.plate_rx, 1)
    if iface.search_rev:
        parse_single(None, gui, iface.reverse_rx, 3)
    iface.rx_fired = True, iface.rx_fired[1]
    reset_sort_size(gui)
    LOG.debug('parse_triple done')


def cell_edit(cell, path, new_text, tree_view, col):
    tree_view.get_model()[path][col] = new_text


def rx_toggle(widget, gui):
    data, iface = gui.data, gui.interface
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
        parse_triple(widget, gui)


def rev_adjust(widget, gui):
    data, iface = gui.data, gui.interface
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
            Gtk.TreeViewColumn(title='reverse', cell_renderer=Gtk.CellRendererText(), text=3), n_cols - 1)
        # cause parsing
        parse_single(None, gui, iface.reverse_rx, 3)
    else:
        iface.search_rev = False
        iface.reverse_rx.set_sensitive(False)
        for col in iface.view_trace_regex.get_columns():
            if col.get_title() == 'reverse':
                iface.view_trace_regex.remove_column(col)

    reset_sort_size(gui)


def reset_sort_size(gui):
    data, iface = gui.data, gui.interface
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        tree_view.columns_autosize()
        for col_index in range(tree_view.get_n_columns()):
            tree_view.get_column(col_index).set_sort_column_id(col_index)

    # check the dataset for empty strings TODO check if stringent enough
    if '' not in commons.get_column(data.trace_store, 2) \
            + commons.get_column(data.trace_store, 4) \
            + commons.get_column(data.plate_store, 0):
        if not iface.plates or '' not in commons.get_column(data.trace_store, 3):
            commons.set_errors(iface, PAGE, False)
            LOG.debug('found no errors')
            return
    assert commons.get_errors(iface, PAGE) or not iface.rx_fired[0] or not iface.rx_fired[1]


def try_online(widget, gui):
    data, iface = gui.data, gui.interface
    LOG.debug('generating online help')
    url = 'https://regex101.com/api/regex/'
    headers = {'content-type': 'application/json'}

    payload = dict(flavor='python', flags='gm', delimiter='"')
    payload['regex'] = iface.single_rx.get_text()
    content = tuple(i.get_text() for i in [iface.wp_rx, iface.single_rx, iface.well_rx,
                                           iface.gene_rx, iface.plate_rx, iface.reverse_rx]) \
              + ('\n'.join([row[1] for row in data.plate_store]),
                 '\n'.join([row[1] for row in data.trace_store]))

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
