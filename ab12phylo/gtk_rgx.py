# 2021 Leo Kaindl

import json
import logging
import random
import re
import string
import threading
import webbrowser
from time import sleep

import gi
import pandas as pd
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo import repo
from ab12phylo_cmd.filter import new_id
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 1
ERRORS = ['groups?', 'no match', 'use groups!']


class rgx_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.single_rt.connect('toggled', self._switch_single_or_triple_regex, )
        iface.triple_rt.connect('toggled', self._switch_single_or_triple_regex)
        iface.triple_rt.join_group(iface.single_rt)
        iface.single_rt.set_active(True)

        iface.regex_apply.connect('clicked', self.parse_all, True)
        iface.wp_apply.connect('clicked', self.parse_single_group, iface.wp_rx, 2)
        iface.single_rx.connect('activate', self.parse_all, True)

        for thing, col in zip(['well', 'gene', 'plate', 'wp'], [2, 4, 3, 2]):
            entry = iface.__getattribute__('%s_rx' % thing)
            entry.connect('activate', self.parse_single_group, entry, col)
            entry.connect('focus_out_event', self.parse_single_group, entry, col)

        iface.view_trace_regex.set_model(data.trace_store)
        iface.view_csv_regex.set_model(data.plate_store)
        iface.reverse_rx_chk.set_active(data.search_rev)
        iface.reverse_rx_chk.connect('toggled', self._switch_search_reverse_reads)
        iface.reverse_rx.connect('activate', self.parse_single_group, iface.reverse_rx, 6)
        iface.regex_try_online.connect('clicked', self.try_online)

        # allow deletion
        for widget, tv in zip([iface.remove_path_regex, iface.remove_csv_regex],
                              [iface.view_trace_regex, iface.view_csv_regex]):
            sel = tv.get_selection()
            sel.set_mode(Gtk.SelectionMode.MULTIPLE)
            widget.connect('clicked', self.delete_files_from_input_selection, PAGE, sel)
            tv.connect('key-press-event', self.tv_keypress, PAGE, sel)

        self.reset_columns()

    def reset_columns(self, do_parse=False):
        data = self.data
        iface = self.iface
        self.check_for_plates()

        # remove old columns:
        for tv in [iface.view_trace_regex, iface.view_csv_regex]:
            [tv.remove_column(col) for col in tv.get_columns()]

        # columns depend on if we are reading plates
        if iface.plates:
            for title, j in zip(['well', 'plate', 'gene', 'file'], [2, 3, 4, 1]):
                crt = Gtk.CellRendererText(editable=True)
                crt.connect('edited', self.cell_edit, iface.view_trace_regex, j)
                col = Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=j, foreground=7)
                col.set_sort_column_id(j)
                iface.view_trace_regex.append_column(col)
            # wellsplates:
            for title, j in zip(['plate ID', 'file'], [2, 1]):
                crm = Gtk.CellRendererText(editable=True)
                crm.connect('edited', self.cell_edit, iface.view_csv_regex, j)
                col = Gtk.TreeViewColumn(title=title, cell_renderer=crm, text=j, foreground=3)
                col.set_sort_column_id(j)
                iface.view_csv_regex.append_column(col)
        else:
            for title, j in zip(['sample', 'gene', 'file'], [2, 4, 1]):
                crt = Gtk.CellRendererText(editable=True)
                crt.connect('edited', self.cell_edit, iface.view_trace_regex, j)
                col = Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=j, foreground=7)
                col.set_sort_column_id(j)
                iface.view_trace_regex.append_column(col)
            # set regex for plate_id column as already fired
            data.rx_fired[-2] = True

        self._switch_search_reverse_reads(self.iface.reverse_rx_chk)
        data.rx_fired[6] = True

        if do_parse:  # fire the initial parse
            self.parse_all(None, force=True)
            self.refresh(page=1)

    def parse_single_group(self, widget, entry, col, fourth=None, refresh_after=True):
        if fourth:  # caught a focus_out_event, therefore shifted arguments
            entry, col = col, fourth
        data = self.data
        iface = self.iface
        LOG.debug('single')
        try:
            regex = re.compile(entry.get_text())
        except re.error:
            self.show_notification('RegEx could not be compiled.')
            return
        errors = self.get_errors(PAGE)
        changed = self.get_changed(PAGE)

        # check if traces or plates
        if widget in [iface.wp_rx, iface.wp_apply]:
            model = data.plate_store
            traces = False
            data.rx_fired[col + data.trace_store.get_n_columns()] = True
        else:
            model = data.trace_store
            traces = True
            data.rx_fired[col] = True

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
                    model[i][-1] = None  # iface.FG
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
                # reverse reads, save boolean values
                try:
                    m = regex.search(file).groups()[0]
                    # reverse read
                    if not changed:
                        changed = row[col] is not True
                    model[i][col] = True
                    if model[i][-1] != iface.RED:
                        model[i][-1] = None  # iface.FG
                except AttributeError as ae:
                    # forward read
                    if not changed:
                        changed = row[col] is not False
                    model[i][col] = False
                    if model[i][-1] != iface.RED:
                        model[i][-1] = None  # iface.FG
                except (IndexError, ValueError):
                    errors, changed = True, True
                    model[i][-1] = iface.RED
                except Exception as ex:
                    raise RuntimeWarning from ex
        self.set_changed(PAGE, changed)
        self.set_errors(PAGE, errors)
        if errors and sum(data.rx_fired) >= 5 and refresh_after:
            self.refresh()

    def parse_all(self, wi, force=False):
        data = self.data
        iface = self.iface
        if not force:  # it's easier to pass an additional param than blocking handlers
            return

        if iface.single_rt.get_active():
            LOG.debug('triple')
            try:
                regex = re.compile(iface.single_rx.get_text())
            except re.error:
                self.show_notification('RegEx could not be compiled')
                return
            errors = self.get_errors(PAGE)
            changed = self.get_changed(PAGE)

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
                        changed = not bool(row[2:5] == [well, plate, gene])
                    data.trace_store[idx] = row[:2] + [well, plate, gene] + row[5:]
                    data.trace_store[idx][-1] = None  # iface.FG
                    continue
                except ValueError as ve:
                    data.trace_store[idx][2:5] = [ERRORS[0], '', '']
                except AttributeError as ae:
                    data.trace_store[idx][2:5] = [ERRORS[1], '', '']
                except IndexError as ie:
                    data.trace_store[idx][2:5] = [ERRORS[2], '', '']
                data.trace_store[idx][7] = iface.RED
                errors, changed = True, True
            self.set_changed(PAGE, changed)
            self.set_errors(PAGE, errors)
            data.rx_fired[2:5] = [True] * 3
        else:
            self.parse_single_group(None, iface.well_rx, 2, False)
            self.parse_single_group(None, iface.gene_rx, 4, False)
        if iface.plates:
            self.parse_single_group(iface.wp_apply, iface.wp_rx, 2, False)
        if iface.reverse_rx_chk.get_active():
            self.parse_single_group(None, iface.reverse_rx, 6, False)
        if self.get_errors(PAGE) and sum(data.rx_fired) >= 5:
            self.refresh()

    def cell_edit(self, cell, path, new_text, tv, col):
        self.save_row_edits(cell, path, new_text, tv, col)
        self.set_changed(PAGE)
        # set row to error-free. also for references
        if tv == self.iface.view_trace_regex:
            mo = tv.get_model()
            mo[path][-1] = self.iface.BLUE if mo[path][5] else None  # iface.FG
        self.refresh()

    def _switch_single_or_triple_regex(self, widget):
        # (In)activates the Entry fields and causes a parse event
        data = self.data
        iface = self.iface
        if widget.get_active():
            self.set_errors(PAGE, False)  # important, do not pre-suppose it will fail
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
            self.parse_all(widget)

    def _switch_search_reverse_reads(self, widget):
        data = self.data
        iface = self.iface
        n_cols = iface.view_trace_regex.get_n_columns()
        if widget.get_active():
            # save in data
            data.search_rev = True
            # enable entry field
            iface.reverse_rx.set_sensitive(True)
            # create new column
            crt = Gtk.CellRendererToggle(radio=False)
            crt.connect('toggled', self._manual_mark_as_reversed)
            col = Gtk.TreeViewColumn(title='rev', cell_renderer=crt, active=6)
            col.set_sort_column_id(6)
            iface.view_trace_regex.insert_column(col, n_cols - 1)
            # set regex as not yet fired
            data.rx_fired[6] = False
            # cause parsing
            self.parse_single_group(None, iface.reverse_rx, 6)
        else:
            # save in data
            data.search_rev = False
            iface.reverse_rx.set_sensitive(False)
            for col in iface.view_trace_regex.get_columns():
                if col.get_title() == 'rev':
                    iface.view_trace_regex.remove_column(col)
            # set regex as already fired
            data.rx_fired[6] = True

    def _manual_mark_as_reversed(self, widget, path):
        """Mark a record as reverse or not"""
        data = self.data
        data.trace_store.set(data.trace_store.get_iter(path), [6], [not data.trace_store[path][6]])
        self.set_changed(PAGE)
        self.refresh()

    def check_for_plates(self):
        iface = self.iface
        plates_now = len(self.data.plate_store) > 0
        if 'plates' not in iface or plates_now != iface.plates:
            iface.plates = plates_now

            # close and inactivate the expander
            iface.wp_expander.set_expanded(iface.plates)
            iface.wp_expander.set_sensitive(iface.plates)
            self.reset_columns()

        if not iface.plates:
            # inactivate the line for the matching single regex
            if not iface.triple_rt.get_active():
                [iface.__getattribute__(name).set_sensitive(False) for name in
                 ['plate_rx', 'plate_regex_label']]

    def refresh(self, **kwargs):
        data = self.data
        iface = self.iface
        self.check_for_plates()
        # make all columns sortable
        for tv in [iface.view_trace_regex, iface.view_csv_regex]:
            tv.columns_autosize()

        # try matching reference files to genes
        data.genes = sorted({gene for gene, is_ref, color in data.trace_store.get_column((4, 5, 7))
                             if not is_ref and color is not iface.RED and gene != ''})
        if not data.genes or data.genes == ['']:
            iface.gene_list.set_label('no genes')
            iface.gene_list.set_visible(False)
            return
        tx = ', '.join(data.genes)
        LOG.debug(f'genes: {tx}')
        iface.gene_list.set_label(tx)
        iface.gene_list.set_visible(True)
        single_gene = data.genes[0] if len(data.genes) == 1 else False

        for i, row in enumerate(data.trace_store):
            # iterate over reference sequences and mark them blue if they are ok
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

        # check the dataset for red lines
        trace_colors = data.trace_store.get_column(-1)
        if iface.RED not in trace_colors and iface.AQUA not in trace_colors \
                and iface.RED not in data.plate_store.get_column(-1):
            self.set_errors(PAGE, False)
            data.rx_fired[2:5] = [True] * 3
            data.rx_fired[6] = True
        else:
            assert self.get_errors(PAGE) or sum(data.rx_fired) < 5

        # check the dataset for versionized entries
        appears = dict()
        for w, p, g in data.trace_store.get_column((2, 3, 4)):
            if w not in appears:
                appears[w] = {p: {g: 1}}
            elif p not in appears[w]:
                appears[w][p] = {g: 1}
            else:
                appears[w][p][g] = appears[w][p].get(g, 0) + 1

        for row in data.trace_store:
            if appears[row[2]][row[3]][row[4]] > 1:
                row[-1] = iface.PURPLE

    def reload_ui_state(self):
        for w_name in ['single_rx', 'well_rx', 'gene_rx', 'plate_rx', 'reverse_rx', 'wp_rx']:
            self.iface.__getattribute__(w_name).set_text(str(self.data.rgx.__getattribute__(w_name)))
        for w_name in ['single_rt', 'triple_rt', 'reverse_rx_chk']:
            self.iface.__getattribute__(w_name).set_active(self.data.rgx.__getattribute__(w_name))

    def try_online(self, *args):
        data = self.data
        iface = self.iface
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
                                'trace filenames:\n%s\n' % content

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
            self.show_notification('online help failed')

    def start_read(self, run_after=None):
        """Called from REFRESH button and shared.proceed"""
        data = self.data
        iface = self.iface

        # stop if there are no traces
        if len(data.trace_store) == 0:
            self.show_notification('No sequence data!')
            return
        if run_after is None:
            self.refresh()  # try matching reference files to genes

        # save_ui_state
        for w_name in ['single_rx', 'well_rx', 'gene_rx', 'plate_rx', 'reverse_rx', 'wp_rx']:
            data.rgx.__setattr__(w_name, iface.__getattribute__(w_name).get_text())
        for w_name in ['single_rt', 'triple_rt', 'reverse_rx_chk']:
            data.rgx.__setattr__(w_name, iface.__getattribute__(w_name).get_active())

        iface.thread = threading.Thread(target=self.do_read)
        iface.run_after = run_after
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        # return to main loop

    def do_read(self):
        data = self.data
        iface = self.iface
        iface.text = ''
        iface.i = 0
        iface.k = len(data.trace_store) + len(data.plate_store) + 2
        data.csvs.clear()
        data.seqdata.clear()
        data.record_order.clear()
        records, warnings, errors = list(), set(), list()

        # use random but reproducible prefixes for references. different seed from raxml -> less confusion
        random.seed(data.seed + 2)

        # dict of reference source organism -> random key
        lookup = dict()

        # read in wellsplates
        LOG.debug('reading wellsplates')
        iface.text = 'reading plates ...'
        for row in data.plate_store:
            df = pd.read_csv(row[0], header=None, engine='python')
            df.index = list(range(1, df.shape[0] + 1))
            df.columns = list(string.ascii_uppercase[0:df.shape[1]])
            file, box = row[1:3]
            if box in data.csvs:
                errors.append('overwrite wellsplate %s with %s' % (box, file))
            data.csvs[box] = df.to_dict()
            iface.i += 1

        # read in FASTAs containing more than one record first
        fasta_paths = {file_path.split('~~')[0] for file_path
                       in data.trace_store.get_column(0) if '~~' in file_path}
        fasta_records = {file_path: {r.id: r for r in SeqIO.parse(
            file_path, 'fasta')} for file_path in fasta_paths}

        # read in trace files
        LOG.debug('reading traces')
        found_numerical = False
        for row in data.trace_store:
            file_path, file, coords, box, gene, is_ref, is_rev, color = row
            iface.text = 'reading %s' % file

            # read in one file
            if file_path.endswith('.ab1') and '~~' not in file_path:
                try:
                    records = [SeqIO.read(file_path, 'abi')]  # ABI traces also only contain a single record!
                except UnicodeDecodeError:
                    errors.append('ABI trace error %s' % file)
            elif '~~' in file_path:
                file_part, record_name = file_path.split('~~')
                records = [fasta_records[file_part][record_name]]
            else:
                try:
                    records = SeqIO.parse(file_path, 'fasta')
                except UnicodeDecodeError:
                    errors.append('Seq file error %s' % file)

            # any kind of seq can define a new gene
            if gene not in data.seqdata:
                data.seqdata[gene] = dict()
                data.metadata[gene] = dict()
                data.genes = sorted(set(data.genes) | {gene})
            keys = list(data.seqdata[gene].keys()) if gene in data.seqdata else list()

            # iterate over records found in file (usually only one)
            for record in records:
                if not is_ref:  # normal sequence
                    try:
                        # swap out well coordinates for isolate numbers
                        (y, x) = (int(coords[1:]), coords[0])
                        _id = data.csvs[box][x][y]
                    except (KeyError, ValueError):
                        if len(data.csvs) == 0:
                            # record.id = record.name.replace(gene, '')
                            if box in ['', ' ', '-', '_']:
                                _id = coords.upper()
                            else:
                                _id = box + '_' + coords
                        else:
                            _id = record.id
                            errors.append('missing %s on wellsplate %s' % (coords, box))
                    attributes = {'file': file_path, 'box': box, 'is_rev': is_rev}

                    if _id.isdigit():
                        # prepend a T to numerical IDs so toytree won't get confused
                        _id = f'T{_id}'
                        found_numerical = True
                    if _id in keys:
                        # add suffix to duplicate IDs
                        warnings.add('duplicate ID %s' % _id)
                        _id = new_id(_id, keys)

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

                    # catch some illegal weird html characters
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
                    attributes = {'file': file_path, 'accession': accession, 'is_rev': is_rev,
                                  'reference_species': species if strain == species
                                  else species + ' strain ' + strain}

                if is_rev:  # strange if it's a reference, but possible anyway
                    record = record.reverse_complement(record.id)
                record = SeqRecord(record.seq, _id, description='', name='',
                                   letter_annotations=getattr(record, 'letter_annotations'))
                # save SeqRecord
                data.seqdata[gene][_id] = record
                data.metadata[gene][_id] = attributes
                # save the order of all records
                data.record_order.append([_id, gene])

            iface.i += 1

        if found_numerical:
            warnings.add('Found records with purely numerical IDs; prepended a T to those.')
        iface.text = 'search missing samples'
        LOG.debug(iface.text)
        data.gene_ids = {g: set(gd.keys()) for g, gd in data.seqdata.items() if g in data.genes}
        all_ids = set.union(*data.gene_ids.values())
        missing = {g: ', '.join(sorted(all_ids - ids))
                   for g, ids in data.gene_ids.items() if all_ids - ids}
        if missing:
            with open(self.wd / repo.PATHS.missing, 'w') as fh:
                fh.write('gene\tmissing samples\n')
                for k, v in missing.items():
                    fh.write('%s\t%s\n' % (k, v))
                    warnings.add('%s missing %s' % (k, v))

        iface.i += 1
        iface.text = 'writing metadata'
        LOG.debug(iface.text)
        self.write_metadata()

        iface.text = 'idle'
        iface.frac = 1
        LOG.debug('reading done')
        GObject.idle_add(self.stop_read, errors, warnings)
        return

    def stop_read(self, errors, warnings):
        """
        Join a thread, then call one or more functions *afterwards*.
        Usually with the same args: gui. These params are stored in iface.
        :param errors:
        :param warnings:
        :return:
        """
        self.iface.thread.join()
        LOG.info('rgx thread idle')
        if errors or warnings:
            self.show_notification('File troubles', ['error: %s' % f for f in errors] +
                                   ['warning: %s' % f for f in warnings], secs=0)
        self.set_changed(PAGE, False)

        # *Now* go back to where you came from
        if self.iface.run_after:
            [do_func() for do_func in self.iface.run_after]
        return
