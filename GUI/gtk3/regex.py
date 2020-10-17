import json
import re
import webbrowser

import gi
import logging
from pathlib import Path
from time import sleep

import requests

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk, Gdk as gdk

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


def init(data, iface):
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
        sel().set_mode(gtk.SelectionMode.MULTIPLE)
        widget.connect('clicked', commons.delete_rows, iface, data, PAGE, sel())

    # connect buttons
    iface.regex_next.connect('clicked', commons.proceed, data, iface)
    iface.regex_back.connect('clicked', commons.step_back, iface)
    reset(data, iface)
    commons.refresh_files(iface, data, PAGE)


def reset(data, iface):
    # create TreeView model
    # path extraction: well/id, (plate), gene, (reversed), file
    data.rx_model = gtk.ListStore(str, str, str, str, str)
    # fill last column with initial filename data
    [data.rx_model.append([''] * 4 + [Path(path[0]).name]) for path in data.trace_model]
    if iface.plates:
        data.wp_model = gtk.ListStore(str, str)
        [data.wp_model.append([''] + [Path(path[0]).name]) for path in data.csv_model]
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
                gtk.TreeViewColumn(title=title, cell_renderer=gtk.CellRendererText(), markup=column))
        # wellsplates:
        for title, column in zip(['plate ID', 'file'], list(range(2))):
            iface.view_csv_regex.append_column(
                gtk.TreeViewColumn(title=title, cell_renderer=gtk.CellRendererText(), markup=column))
    else:
        for title, column in zip(['sample', 'gene', 'file'], [0, 2, 4]):
            iface.view_trace_regex.append_column(
                gtk.TreeViewColumn(title=title, cell_renderer=gtk.CellRendererText(), markup=column))

    reset_sort_size(iface)
    iface.reverse_rx_chk.set_active(False)


def parse_single(widget, data, iface, entry, col):
    LOG.debug('parsing from %s' % entry.get_name())
    regex = re.compile(entry.get_text())

    if widget in [iface.wp_rx, iface.wp_apply]:
        model = data.wp_model
    else:
        model = data.rx_model

    if entry is not iface.reverse_rx:
        for row in range(len(model)):
            file = model[row][-1]
            try:
                model[row][col] = regex.search(file).groups()[0]
            except ValueError as ve:
                # maybe wrong number of groups
                model[row][col] = ERRORS[0]
            except AttributeError as ae:
                # no match
                model[row][col] = ERRORS[1]
            except IndexError as ie:
                # no groups used
                model[row][col] = ERRORS[2]
    else:
        for row in range(len(model)):
            file = model[row][-1]
            try:
                a = regex.search(file).groups()[0]
                model[row][col] = MARKUP[0]
            except ValueError as ve:
                model[row][col] = MARKUP[1]
            except AttributeError as ae:
                model[row][col] = MARKUP[2]
            except IndexError as ie:
                model[row][col] = MARKUP[3]
            # try:
            #     a = regex.search(file).groups()[0]
            #     model[row][col] = 'go-previous-symbolic'
            # except ValueError:
            #     model[row][col] = 'gtk-dialog-warning'
            # except  AttributeError:
            #     model[row][col] = ''  # 'go-next-symbolic'
            # except IndexError:
            #     model[row][col] = 'gtk-dialog-warning'


def parse_triple(widget, data, iface):
    if iface.single_rt.get_active():
        LOG.debug('parsing with single regex')
        # if parsing with only a single regex
        regex = re.compile(iface.single_rx.get_text())

        for row in range(len(data.rx_model)):
            file = data.rx_model[row][-1]
            try:
                m = regex.search(file)
                plate, gene, well = m.groups() if iface.plates else (None, *m.groups())
                data.rx_model[row] = [well, plate, gene, data.rx_model[row][3], file]
            except ValueError as ve:
                data.rx_model[row][0] = ERRORS[0]
            except AttributeError as ae:
                data.rx_model[row][0] = ERRORS[1]
            except IndexError as ie:
                data.rx_model[row][0] = ERRORS[2]
    else:
        parse_single(None, data, iface, iface.well_rx, 0)
        parse_single(None, data, iface, iface.gene_rx, 2)
        if iface.plates:
            parse_single(None, data, iface, iface.plate_rx, 1)
    if iface.search_rev:
        parse_single(None, data, iface, iface.reverse_rx, 3)
    reset_sort_size(iface)
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
        # px = gtk.CellRendererPixbuf()
        # col = gtk.TreeViewColumn(title=' ', cell_renderer=px)
        # col.set_cell_data_func(px,
        #                        lambda column, cell, model, _iter, user_data:
        #                        cell.set_property('icon-name', model[_iter][3]))
        # iface.view_trace_regex.insert_column(col, iface.view_trace_regex.get_n_columns() - 1)
        iface.view_trace_regex.insert_column(
            gtk.TreeViewColumn(title='reverse', cell_renderer=gtk.CellRendererText(), markup=3), n_cols - 1)
        # cause parsing
        parse_single(None, data, iface, iface.reverse_rx, 3)
    else:
        iface.search_rev = False
        iface.reverse_rx.set_sensitive(False)
        for col in iface.view_trace_regex.get_columns():
            if col.get_title() == 'reverse':
                iface.view_trace_regex.remove_column(col)

    reset_sort_size(iface)


def reset_sort_size(iface):
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        tree_view.columns_autosize()
        for col_index in range(tree_view.get_n_columns()):
            tree_view.get_column(col_index).set_sort_column_id(col_index)


def try_online(widget, data, iface):
    LOG.debug('generating online help')
    url = 'https://regex101.com/api/regex/'
    headers = {'content-type': 'application/json'}

    payload = dict(flavor='python', flags='gm', delimiter='"')
    payload['regex'] = iface.single_rx.get_text()
    content = tuple(i.get_text() for i in [iface.wp_rx, iface.single_rx, iface.well_rx,
                                           iface.gene_rx, iface.plate_rx, iface.reverse_rx]) \
              + ('\n'.join([row[-1] for row in data.wp_model]),
                 '\n'.join([row[-1] for row in data.rx_model]))

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
