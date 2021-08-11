# 2021 Leo Kaindl

import logging
import shutil
import stat
import subprocess
import threading
from argparse import Namespace
from pathlib import Path
from time import sleep

import gi
import pandas as pd
from Bio import SeqIO
from numpy import isnan, nan

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo_cmd import blast
from ab12phylo import repo
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 5


class blast_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.blast_spinner.stop()
        iface.spinbox.set_visible(False)

        iface.remote_db.set_model(data.remote_dbs)
        iface.remote_db.get_child().connect('focus_out_event', self._save_custom_remote_db)

        # connect buttons
        iface.blast_local.connect('clicked', self.start_BLAST, 'local')
        iface.blast_remote.connect('clicked', self.start_BLAST, 'remote')
        iface.blast_import.connect('clicked', self.start_BLAST, None)
        iface.xml_import.connect('file_set', self.start_BLAST, None)
        iface.blast_exe.connect('file_set', lambda *args: self.prep0(iface.blast_exe.get_filename()))
        iface.blast_db.connect('button-release-event', lambda co, *args: co.
                               popdown() if co.props.popup_shown else co.popup())
        iface.blast_gene.connect('changed', lambda *args: self._refill())

        # set up the species annotation table
        spi = iface.sp_info
        spi.connect_after('size-allocate', lambda *args: [col.
                          set_max_width(spi.get_allocated_width() // 2) for col in spi.get_columns()])
        spi.set_model(data.sp_model)
        spi.set_reorderable(True)
        for i, title in enumerate(['sample', 'pid', 'species', 'runner-up']):
            crt = Gtk.CellRendererText(editable=True)
            crt.connect('edited', self._save_sp_edit, spi, i)
            col = Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=i, foreground=4)
            col.set_sort_column_id(i)
            col.set_resizable(True)
            spi.append_column(col)
        spi.columns_autosize()

        iface.blast_seen = False

    def prep0(self, path):
        data = self.data
        iface = self.iface
        iface.i = 0
        iface.k = 3
        iface.text = 'search for BLAST+ installation'
        if path:
            path = Path(path).resolve()
            if path.is_file():
                path = path.parent  # path is supposed to be a directory!
            data.blast_path = path  # save path to executable in project

        # check if all necessary BLAST+ scripts are available
        try:
            binary = ab12phylo_app_base.CFG['blastn']
        except KeyError:
            binary = shutil.which('blastn', path=path)
        if binary is None:
            binary = shutil.which('blastn')

        scripts = ['blastn', 'blastdbcmd', 'update_blastdb.pl', 'makeblastdb']
        if not binary:
            missing = scripts
        else:
            missing = list()
            binary = Path(binary)

            for sc in scripts:
                try:
                    exe = next(binary.parent.rglob(f'{sc}{repo.EXE if not sc.endswith(".pl") else ""}'))
                    # Make the file executable
                    try:
                        exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
                    except:
                        pass
                except StopIteration:
                    missing.append(sc)

        if not missing:
            data.blast_path = binary.parent
            ab12phylo_app_base.CFG['blastn'] = binary
        else:
            data.blast_path = None
            self.show_notification(msg='Necessary scripts could not be '
                                       'found on the $PATH:', items=missing)
            sleep(.1)
            GObject.idle_add(self.prep4)
            return True

        iface.blast_exe.set_filename(str(data.blast_path))
        mo = Gtk.ListStore(str, str, str, str, bool)
        iface.blast_db.set_model(mo)  # shared model
        iface.db_info.set_model(mo)
        iface.db_info.set_reorderable(True)
        for i, title in enumerate(['database', 'size', 'title']):
            col = Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(),
                                     text=i + 1, underline=4)
            col.set_resizable(True)
            iface.db_info.append_column(col)

        iface.thread = threading.Thread(target=self.prep1, args=[data.blast_path, mo])
        GObject.timeout_add(50, self.update, PAGE)
        iface.i = 1
        iface.thread.start()
        return True

    def prep1(self, path, mo):
        """
        There is a BLAST+ executable on the $PATH. Pre-select it, search for local and remote
        databases with BLAST+ commands, then set-up the db_info and the blast_db ComboBox.
        :param gui:
        :param binary:
        :return:
        """
        data = self.data
        iface = self.iface
        try:
            iface.text = 'look for local databases'
            output = subprocess.check_output(
                str(path / 'blastdbcmd') + ' -recursive -list_outfmt "%f\t%t\t%U" -list ' +
                str(repo.BASE_DIR), shell=True).decode('utf-8').strip().split('\n')
            if output and output != ['']:
                sleep(.05)
                GObject.idle_add(self.prep2, iface.db_info, output, mo)
                iface.i = 2
        except subprocess.CalledProcessError as ex:
            LOG.error(ex)

        try:
            iface.text = 'look for online pre-compiled databases'
            update_arg = f"{shutil.which('perl')} {path / 'update_blastdb.pl'} --source gcp --showall tsv"
            LOG.debug(update_arg)
            output = subprocess.check_output(update_arg, shell=True).decode('utf-8')
            sleep(.05)
            GObject.idle_add(self.prep3, data.remote_dbs, output, mo, iface.remote_db)
            iface.i = 3
        except subprocess.CalledProcessError as ex:
            LOG.error(ex)
        sleep(.05)
        GObject.idle_add(self.prep4)
        return True

    @staticmethod
    def prep2(db_info, output, mo):
        db_info.set_tooltip_text('Underlining indicates a local database')
        for db in output:
            db = db.split('\t')
            try:
                mo.append([db[0], Path(db[0]).name, human_bytes(
                    sum(f.stat().st_size for f in Path(db[0]).parent
                        .glob('**/*') if f.is_file())), db[1], True])
            except IndexError:
                pass  # if directory does not exist or is empty

    @staticmethod
    def prep3(remote_dbs, output, mo, remote_db):
        remote_dbs.clear()
        for i, db in enumerate(output.strip().split('\n')[1:]):
            db = db.split('\t')
            # fill db_info table
            mo.append([None, db[0], human_bytes(
                2 ** 30 * float(db[2].strip())), db[1], False])
            # fill remote_db ComboBox
            remote_dbs.append([db[0], i])
            if db[0] == 'nt':
                remote_db.set_active(i)
        if len(remote_dbs) == 0:
            remote_dbs.append(['no remotes found', 0])

    def prep4(self):
        iface = self.iface
        if iface.thread.is_alive():
            iface.thread.join()
        # sync db selection in db_info TreeView and blast_db ComboBox
        sel = iface.db_info.get_selection()

        def connector(selection, combobox):
            try:
                combobox.set_active(selection.get_selected_rows()[1][0].get_indices()[0])
            except IndexError:
                pass

        sel.connect('changed', connector, iface.blast_db)
        sel.select_path(0)
        iface.blast_db.connect('changed', lambda *args: sel
                               .select_path(iface.blast_db.get_active()))
        LOG.info('BLAST prep done')
        return False

    def refresh(self):
        """Re-view the page. Re-load the current genes to the marker-gene selector"""
        data = self.data
        iface = self.iface

        if not iface.blast_seen:
            # look for BLAST+ executable in the config and on the $PATH; fill db_info table
            for path in [shutil.which('blastn'), data.blast_path]:
                if self.prep0(path):
                    break
            iface.blast_seen = True

        iface.blast_gene.remove_all()
        genes = list(data.genes)
        [iface.blast_gene.append_text(gene) for gene in genes]
        if data.gene_for_preview and data.gene_for_preview in genes:
            idx = genes.index(data.gene_for_preview)
        else:
            idx = 0
        iface.blast_gene.set_active(idx)

        # always allow skipping BLAST
        self.set_changed(PAGE, False)

    def start_BLAST(self, widget, mode, *args):
        """Set-up the BLAST thread"""
        data = self.data
        iface = self.iface

        if 'blast_wrapper' in iface:
            # The button was pressed when it was in the 'Stop' state
            iface.pill2kill.set()
            sleep(.1)
            if iface.blaster.is_alive():
                sleep(.1)
                iface.blaster.stop()
            return

        iface.blast_spinner.start()
        iface.spinbox.set_visible(True)
        gene = iface.blast_gene.get_active_text()
        data.gene_for_preview = gene
        LOG.info('running %s for %s' % (mode, gene))
        if mode == 'local':
            fasta = self.wd / gene / (gene + '.fasta')
            if not fasta.is_file():
                self.show_notification('%s with sequences for BLAST does not exist' % fasta)
                iface.blast_spinner.stop()
                iface.spinbox.set_visible(False)
                return
            pars = [True, False, False]
            adapt_button = iface.blast_local, mode
            # re-read all sequence data
            seqdata = {g: {r.id: r for r in SeqIO.parse(
                self.wd / g / (g + '_all.fasta'), 'fasta')} for g in data.genes}

        elif mode == 'remote':
            fasta = iface.missing_fasta.get_filename()
            if not fasta and mode == 'remote':
                self.show_notification('select a FASTA with sequence data first')
                iface.blast_spinner.stop()
                iface.spinbox.set_visible(False)
                return
            pars = [False, True, False]
            adapt_button = iface.blast_remote, mode
            # read seqdata from selected file
            seqdata = {gene: {r.id: r for r in SeqIO.parse(fasta, 'fasta')}}

        else:  # XML import
            pars = [True, True, iface.xml_import.get_filenames()]
            adapt_button = iface.blast_import, 'Load'
            # re-read all sequence data
            seqdata = {g: {r.id: r for r in SeqIO.parse(
                self.wd / g / (g + '_all.fasta'), 'fasta')} for g in data.genes}

        pars = dict(zip(['no_remote', 'no_local', 'BLAST_xml'], pars))
        (im, la), tx = adapt_button[0].get_child().get_children(), adapt_button[1]
        im.set_from_icon_name('media-playback-stop-symbolic', 4)
        la.set_text('Stop')
        db_row = iface.blast_db.get_model()[iface.blast_db.get_active()]

        # define static parameters
        args = {'gui': True,
                'cfg': ab12phylo_app_base.CFG,
                'no_BLAST': False,
                'genes': [gene],
                'db': db_row[1],
                'remote_db': data.remote_dbs[iface.remote_db.get_active()][0],
                'timeout': repo.DOWNLOAD_TIMEOUT,
                'dir': self.wd,
                'xml': self.wd / repo.PATHS.xml,
                'www_xml': str(self.wd / repo.PATHS.www_xml),
                'tsv': self.wd / repo.PATHS.tsv,
                'bad_seqs': self.wd / repo.PATHS.bad_seqs,
                'missing_fasta': self.wd / repo.PATHS.missing_fasta,
                'dbpath': Path(db_row[0]).parent if db_row[0] else None}
        args.update(pars)
        args = Namespace(**args)
        reader = Namespace(**{'seqdata': seqdata, 'metadata': None})

        iface.blast_tup = (im, la, tx, gene)
        LOG.debug('init blast.blast_build with %s' % str(args))
        iface.blaster = blast.blast_build(args, reader)
        iface.blast_wrapper = threading.Thread(
            # iface.blaster = multiprocessing.Process(
            target=self.do_BLAST)
        GObject.timeout_add(1000, self.update_BLAST)
        iface.blast_wrapper.start()
        return
        # return to main loop

    def do_BLAST(self):
        """Run BLAST thread"""
        try:
            self.iface.blaster.run()
        except Exception as ex:
            self.show_notification(f'BLAST failed: {ex}')
        GObject.idle_add(self.stop_BLAST)

    def stop_BLAST(self):
        """Finish the BLAST thread"""
        data = self.data
        iface = self.iface
        if 'blaster' in iface:
            blaster = iface.blaster
            if blaster.update:  # a local database was created
                iface.blast_db.get_model()[iface.blast_db.get_active()][-1] = True
            if iface.blaster.missing_fasta.is_file():
                iface.missing_fasta.set_filename(str(blaster.missing_fasta))
            if blaster.TSV != self.wd / repo.PATHS.tsv:
                shutil.move(blaster.TSV, self.wd / repo.PATHS.tsv)
                # df_old = self._df_with_sp(iface.blaster.TSV)
                # df_old = df_old.loc[df_old['BLAST_species'] != '', :]
                # df_new = self._df_with_sp(self.wd / repo.PATHS.tsv)
                # for _id, series in df_old.iterrows():
                #     df_new.loc[(df_new.index == series.name) & (df_new['gene'] == series.gene), :] = \
                #         df_old.loc[(df_old.index == series.name) & (df_old['gene'] == series.gene), :]

        im, la, tx, gene = iface.blast_tup
        self._refill(gene, fresh=True)

        iface.blast_spinner.stop()
        iface.spinbox.set_visible(False)
        im.set_from_icon_name('media-playback-start-symbolic', 4)
        la.set_text(tx)

        if iface.pill2kill.is_set():
            tx = 'stopped BLAST'
            self.show_notification(msg=tx, secs=2)
        else:
            tx = 'BLAST finished'
            iface.blast_help_stack.set_visible_child_name('sp_info')
            if iface.notebook.get_current_page() != PAGE:
                self.show_notification(msg=tx, secs=5)
            else:
                self.save(silent=True)
        LOG.info(tx)
        iface.pill2kill.clear()
        del iface.blast_wrapper, iface.blaster
        self.set_changed(PAGE, False)

    @staticmethod
    def _save_custom_remote_db(entry, *args):
        combo = entry.get_parent().get_parent()
        if not combo.get_active_iter():
            tx = entry.get_text()
            combo.get_model().append([tx, -1])
            LOG.debug('entered remote_db %s' % tx)

    @staticmethod
    def _df_with_sp(path):
        LOG.debug('reading df')
        df = pd.read_csv(path, sep='\t', dtype={'id': str})
        df = df.set_index('id')
        if 'BLAST_species' in df:
            df.BLAST_species.fillna('', inplace=True)
        if 'extra_species' in df:
            df.extra_species.fillna('', inplace=True)
        if 'pid' not in df:
            df['pid'] = nan
        return df

    def _refill(self, gene=None, fresh=False):
        """Re-fill the species annotation table"""
        data = self.data
        iface = self.iface
        if not gene:
            gene = iface.blast_gene.get_active_text()
        LOG.debug('re-fill species annotations for %s' % gene)
        if 'df' not in iface.tempspace or fresh:
            try:
                iface.tempspace.df = self._df_with_sp(self.wd / repo.PATHS.tsv)
            except FileNotFoundError as ex:
                LOG.debug(ex)
        df = iface.tempspace.df
        data.sp_model.clear()
        df = df.loc[df['gene'] == gene]
        for sample, r in df.iterrows():
            if 'pid' not in r or isnan(r.pid):
                pid, r.BLAST_species, r.extra_species = '', '', ''
            elif type(r.pid) == float:
                pid = '%.2f' % r.pid if r.pid < 100 else '100'
            else:
                pid = r.pid
            color = None
            if 'quality' in r and not pd.isna(r.quality) and r.quality != 'no phreds':
                color = iface.AQUA
            # and again:
            row = [str(i) for i in [sample, pid, r.BLAST_species, r.extra_species]] + [color]
            data.sp_model.append(row)

    def _save_sp_edit(self, cell, path, new_text, tv, col):
        iface = self.iface
        if col == 0:
            return  # do not allow sample IDs to be changed
        if 'df' not in iface.tempspace:
            iface.tempspace.df = self._df_with_sp(self.wd / repo.PATHS.tsv)
        df = iface.tempspace.df
        mo = tv.get_model()
        gene = iface.blast_gene.get_active_text()
        c2 = {1: 'pid', 2: 'BLAST_species', 3: 'extra_species'}
        if col == 1:
            try:
                new_text = float(new_text)
            except ValueError:
                new_text = 0
        # df.loc[df['gene'] == gene].loc[mo[path][0]][c2[col]] = new_text
        # pandas DataFrame get with multiple conditions including index then assign
        df.loc[(df.index == mo[path][0]) & (df['gene'] == gene), c2[col]] = new_text
        self.save_row_edits(cell, path, str(new_text), tv, col)
        self.set_changed(PAGE)


def human_bytes(num):
    step_unit = 1024.0
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < step_unit:
            return '{:.0f} {}'.format(num, x)
        num /= step_unit


class bump_log_level:

    def __init__(self, log, off=False):
        self.level = max(log.level, logging.INFO)
        self.off = off

    def __enter__(self):
        if not self.off:
            logging.disable(self.level)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        logging.disable(logging.NOTSET)
