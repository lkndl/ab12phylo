import json
import logging
import re
import string
import sys
import threading
import webbrowser
from pathlib import Path
from time import sleep

import gi
import pandas as pd
import requests, random
from Bio import SeqIO

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import commons, quality
from ab12phylo import filter

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 1
ERRORS = ['groups?', 'no match', 'use groups!']


def init(gui):
    data, iface = gui.data, gui.interface
    iface.plates = True

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

    # iface.rx_fired stores if these columns have been parsed before:
    # [bool] with trace and then plate columns. sum(iface.rx_fired) = 5
    iface.rx_fired = [False] * (data.trace_store.get_n_columns() + data.plate_store.get_n_columns())

    # columns depend on if we are reading plates
    if iface.plates:
        for title, column in zip(['well', 'plate', 'gene', 'file'], [2, 3, 4, 1]):
            crt = Gtk.CellRendererText(editable=True)
            crt.connect('edited', cell_edit, iface, iface.view_trace_regex, column)
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=column, foreground=7))
        # wellsplates:
        for title, column in zip(['plate ID', 'file'], [2, 1]):
            crm = Gtk.CellRendererText(editable=True)
            crm.connect('edited', cell_edit, iface, iface.view_csv_regex, column)
            iface.view_csv_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crm, text=column, foreground=3))
    else:
        for title, column in zip(['sample', 'gene', 'file'], [2, 4, 1]):
            crt = Gtk.CellRendererText(editable=True)
            crt.connect('edited', cell_edit, iface, iface.view_trace_regex, column)
            iface.view_trace_regex.append_column(
                Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=column, foreground=7))
        # set regex for plate_id column as already fired
        iface.rx_fired[-2] = True

    # reset_sort_size(gui)
    iface.reverse_rx_chk.set_active(False)
    iface.rx_fired[6] = True

    # fire the initial parse
    parse_single(iface.wp_rx, gui, iface.wp_rx, 2)
    parse_triple(None, gui)


def parse_single(widget, gui, entry, col, fifth=None):
    if fifth:  # caught a focus_out_event, therefore shifted arguments
        gui, entry, col = entry, col, fifth

    data, iface = gui.data, gui.interface
    LOG.debug('parsing from %s' % entry.get_name())
    try:
        regex = re.compile(entry.get_text())
    except re.error:
        commons.show_message_dialog('RegEx could not be compiled.')
        return
    errors = commons.get_errors(iface, PAGE)
    changed = commons.get_changed(iface, PAGE)

    # check if traces or plates
    if widget in [iface.wp_rx, iface.wp_apply]:
        model = data.plate_store
        traces = False
        iface.rx_fired[col + data.trace_store.get_n_columns()] = True
    else:
        model = data.trace_store
        traces = True
        iface.rx_fired[col] = True

    # parse each filename
    for i, row in enumerate(model):
        # skip references
        if traces and row[5]:
            errors = errors or row[-1] == iface.AQUA
            continue

        file = row[1]
        if entry is not iface.reverse_rx:
            try:
                m = regex.search(file).groups()[0]
                if not changed:
                    changed = m != row[col]
                model[i][col] = m
                model[i][-1] = iface.FG
                continue
            except ValueError as ve:
                # maybe wrong number of groups
                model[i][col] = ERRORS[0]
            except AttributeError as ae:
                # no match
                model[i][col] = ERRORS[1]
            except IndexError as ie:
                # no groups used
                model[i][col] = ERRORS[2]
            model[i][-1] = iface.RED
            errors, changed = True, True
        else:
            # MARK reverse reads, save boolean values
            try:
                m = regex.search(file).groups()[0]
                # reverse read
                if not changed:
                    changed = row[col] is not True
                model[i][col] = True
                model[i][-1] = iface.FG
            except AttributeError as ae:
                # forward read
                if not changed:
                    changed = row[col] is not False
                model[i][col] = False
                model[i][-1] = iface.FG
            except (IndexError, ValueError):
                errors, changed = True, True
                model[i][-1] = iface.RED
            except Exception:
                assert False
    commons.set_changed(iface, PAGE, changed)
    commons.set_errors(iface, PAGE, errors)
    if entry is iface.gene_rx:
        search_genes(gui)
    re_check(gui)


def parse_triple(widget, gui):
    data, iface = gui.data, gui.interface
    if iface.single_rt.get_active():
        LOG.debug('parsing with single regex')
        try:
            regex = re.compile(iface.single_rx.get_text())
        except re.error:
            commons.show_message_dialog('RegEx could not be compiled')
            return
        errors = commons.get_errors(iface, PAGE)
        changed = commons.get_changed(iface, PAGE)

        for idx, row in enumerate(data.trace_store):
            # skip references
            if row[5]:  # is reference
                errors = errors or row[-1] == iface.AQUA
                continue
            file = row[1]
            try:
                m = regex.search(file)
                plate, gene, well = m.groups() if iface.plates else (None, m.groups())
                if not changed:
                    # TODO continue here
                    changed = not bool(row == row[:2] + [well, plate, gene] + row[5:])
                    print('what\'cha gonna do?', file=sys.stderr)
                data.trace_store[idx] = row[:2] + [well, plate, gene] + row[5:]
                data.trace_store[idx][-1] = iface.FG
                continue
            except ValueError as ve:
                data.trace_store[idx][2:5] = [ERRORS[0], '', '']
            except AttributeError as ae:
                data.trace_store[idx][2:5] = [ERRORS[1], '', '']
            except IndexError as ie:
                data.trace_store[idx][2:5] = [ERRORS[2], '', '']
            data.trace_store[idx][-1] = iface.RED
            errors, changed = True, True
        commons.set_changed(iface, PAGE, changed)
        commons.set_errors(iface, PAGE, errors)
        iface.rx_fired[2:5] = [True] * 3
        search_genes(gui)
        re_check(gui)
    else:
        parse_single(None, gui, iface.well_rx, 2)
        parse_single(None, gui, iface.gene_rx, 4)
        if iface.plates:
            parse_single(None, gui, iface.plate_rx, 3)
    if iface.reverse_rx_chk.get_active():
        parse_single(None, gui, iface.reverse_rx, 6)
    LOG.debug('parse_triple done')


def cell_edit(cell, path, new_text, iface, tv, col):
    mo = tv.get_model()
    old_text = mo[path][col]
    if old_text == new_text:
        return
    mo[path][col] = new_text
    commons.set_changed(iface, PAGE, True)
    # set row to error-free. also for references
    if tv == iface.view_trace_regex:
        mo[path][-1] = iface.BLUE if mo[path][5] else iface.FG


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
        # enable entry field
        iface.reverse_rx.set_sensitive(True)
        # create new column
        iface.view_trace_regex.insert_column(
            Gtk.TreeViewColumn(title='rev', cell_renderer=Gtk.CellRendererToggle(radio=False), active=6), n_cols - 1)
        # set regex as not yet fired
        iface.rx_fired[6] = False
        # cause parsing
        parse_single(None, gui, iface.reverse_rx, 6)
    else:
        iface.reverse_rx.set_sensitive(False)
        for col in iface.view_trace_regex.get_columns():
            if col.get_title() == 'rev':
                iface.view_trace_regex.remove_column(col)
        # set regex as already fired
        iface.rx_fired[6] = True


def search_genes(gui):
    # try matching reference files to genes
    data, iface = gui.data, gui.interface

    data.genes = {gene for gene, is_ref, color in commons.get_column(data.trace_store, (4, 5, 7))
                  if not is_ref and color is not iface.RED}
    if not data.genes or data.genes == {''}:
        return
    single_gene = data.genes.pop() if len(data.genes) == 1 else False
    for i, row in enumerate(data.trace_store):
        if row[5] and row[-1] == iface.AQUA:
            if single_gene:
                data.trace_store[i][4] = single_gene
                data.trace_store[i][7] = iface.BLUE
            else:
                file_name = row[1].upper()
                for gene in data.genes:
                    if gene.upper() in file_name:
                        data.trace_store[i][4] = gene
                        data.trace_store[i][7] = iface.BLUE
                        break


def re_check(gui):
    data, iface = gui.data, gui.interface
    # make all columns sortable
    for tree_view in [iface.view_trace_regex, iface.view_csv_regex]:
        tree_view.columns_autosize()
        for col_index in range(tree_view.get_n_columns()):
            tree_view.get_column(col_index).set_sort_column_id(col_index)

    # check the dataset for red lines
    if iface.RED not in commons.get_column(data.trace_store, -1) \
            and iface.AQUA not in commons.get_column(data.trace_store, - 1) \
            and iface.RED not in commons.get_column(data.plate_store, -1):
        commons.set_errors(iface, PAGE, False)
        LOG.debug('found no errors')
    else:
        LOG.debug('found errors')
        assert commons.get_errors(iface, PAGE) or sum(iface.rx_fired) < 5


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


def read_files(gui):
    data, iface = gui.data, gui.interface

    # stop if there are no traces
    if len(data.trace_store) == 0:
        commons.show_message_dialog('No sequence data!')
        # return  # TODO maybe really stop

    iface.thread = threading.Thread(target=read, args=[gui])
    iface.running = True
    GObject.timeout_add(100, commons.update, iface, iface.read_prog, PAGE)
    iface.thread.start()
    # GUI thread returns to main loop


def read(gui):
    data, iface = gui.data, gui.interface
    iface.frac = 0
    iface.txt = ''
    data.csvs.clear()
    data.seqdata.clear()
    data.record_order.clear()
    records, warnings, errors = list(), set(), set()
    all_there_is_to_do = len(data.trace_store) + len(data.plate_store)
    done = 0

    # use random but reproducible prefixes for references. different seed from raxml -> less confusion
    random.seed(data.seed + 2)

    # dict of reference source organism -> random key
    lookup = dict()

    # read in wellsplates
    LOG.debug('reading wellsplates')
    iface.txt = 'reading plates ...'
    for row in data.plate_store:
        df = pd.read_csv(row[0], header=None, engine='python')
        df.index = list(range(1, df.shape[0] + 1))
        df.columns = list(string.ascii_uppercase[0:df.shape[1]])
        file, box = row[1:3]
        if box in data.csvs:
            errors.add('overwrite wellsplate %s with %s' % (box, file))
        data.csvs[box] = df
        done += 1
        iface.frac = done / all_there_is_to_do

    # read in trace files
    LOG.debug('reading traces')
    for row in data.trace_store:
        file_path, file, coords, box, gene, is_ref, is_rev, color = row
        iface.txt = 'reading %s' % file

        # read in one file
        if file_path.endswith('.ab1'):
            try:
                records = [SeqIO.read(file_path, 'abi')]  # ABI traces also only contain a single record!
            except UnicodeDecodeError:
                errors.add('ABI trace error %s' % file)
        elif file_path.endswith('.fasta') or file_path.endswith('.fa') or file_path.endswith('.seq'):
            try:
                records = SeqIO.parse(file_path, 'fasta')
            except UnicodeDecodeError:
                errors.add('Seq file error %s' % file)

        # iterate over records found in file (usually only one)
        for record in records:

            # any kind of seq can define a new gene
            if gene not in data.seqdata:
                data.seqdata[gene] = dict()
                data.metadata[gene] = dict()
                data.genes.add(gene)
            elif record.id in data.seqdata[gene]:
                warnings.add('duplicate ID %s' % record.id)
                # add suffix to duplicate IDs
                record = filter.new_version(record, data.seqdata[gene].keys())

            # normal sequence
            if not is_ref:
                try:
                    # swap out well coordinates for isolate numbers
                    (y, x) = (int(coords[1:]), coords[0])
                    record.id = data.csvs[box].loc[y, x]
                except (KeyError, ValueError):
                    if len(data.csvs) == 0:
                        # record.id = record.name.replace(gene, '')
                        if box in ['', ' ', '-', '_']:
                            record.id = coords.upper()
                        else:
                            record.id = box + '_' + coords
                    else:
                        errors.add('missing wellsplate %s' % box)

                attributes = {'file': file_path, 'wellsplate': box, 'is_ref': False, 'is_rev': is_rev}
                # TODO continue here

                if is_rev:
                    record = record.reverse_complement(record.id, description='')

            else:  # reference
                # parse species and possibly strain
                strain = re.split(r'[\s_]', record.description.strip().split(',')[0])
                accession = strain.pop(0)
                species = strain.pop(0) + ' ' + strain.pop(0)
                try:
                    ix = strain.index('strain')
                    strain = ' '.join(strain[ix + 1:ix + 3])
                    # cope with longer composite species names
                    if ix > 0:
                        species += ' ' + ' '.join(strain[0:ix])
                except ValueError:
                    strain = species

                # catch some illegal characters
                for char in ['<', '>', '\'', '"', '&']:
                    strain = strain.replace(char, '')

                # retrieve or generate random key
                if strain in lookup:
                    _id = lookup[strain]
                else:
                    # swap original ids for random short ones that no tool takes offense at.
                    _id = 'REF_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
                    lookup[strain] = _id

                # save original id+description
                attributes = {'file': file_path, 'accession': accession, 'is_rev': False, 'is_ref': True,
                              'reference_species': species + ' strain ' + strain}
                record.id = _id
                record.description = ''  # MARK do not delete deletion

            # save SeqRecord
            data.seqdata[gene][record.id] = record
            data.metadata[gene][record.id] = attributes
            # save the order of all records
            data.record_order.append((record.id, gene))

        done += 1
        iface.frac = done / all_there_is_to_do

    LOG.debug('reading done')
    GObject.idle_add(stop, gui, errors, warnings)
    return


def stop(gui, errors, warnings):
    data, iface = gui.data, gui.interface
    iface.running = False
    iface.thread.join()
    iface.read_prog.set_text('idle')
    LOG.info('idle')
    if errors:
        commons.show_message_dialog('There were errors reading some files', errors)
    if warnings:
        commons.show_message_dialog('Additional warnings', warnings)
    # now finally flip to next page
    commons.set_changed(iface, PAGE, False)
    commons.proceed(None, gui)
    quality.reset(gui)
    return
