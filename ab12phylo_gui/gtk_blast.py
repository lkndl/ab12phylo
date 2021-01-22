# 2020 Leo Kaindl

import logging
import shutil
import subprocess
import threading
from argparse import Namespace
from pathlib import Path
from time import sleep

import gi
import pandas as pd
from Bio import SeqIO
from numpy import isnan, nan

from ab12phylo import blast

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo_gui import shared
from ab12phylo_gui.static import PATHS, BASE_DIR, DOWNLOAD_TIMEOUT

LOG = logging.getLogger(__name__)
PAGE = 5


def init(gui):
    """Initialize the page. Connect buttons"""
    data, iface = gui.data, gui.iface

    iface.pill2kill = threading.Event()
    iface.blast_spinner.stop()
    iface.spinbox.set_visible(False)

    iface.remote_db.set_model(data.remote_dbs)
    iface.remote_db.get_child().connect('focus_out_event', _save_custom_remote_db)

    # connect buttons
    iface.blast_local.connect('clicked', start_BLAST, gui, 'local')
    iface.blast_remote.connect('clicked', start_BLAST, gui, 'remote')
    iface.blast_import.connect('clicked', start_BLAST, gui, None)
    iface.xml_import.connect('file_set', start_BLAST, gui, None)
    iface.blast_exe.connect('file_set', lambda *args: prep0(gui, iface.blast_exe.get_filename()))
    iface.blast_db.connect('button-release-event', lambda co, *args: co.
                           popdown() if co.props.popup_shown else co.popup())
    iface.blast_gene.connect('changed', lambda *args: _refill(gui))  # TODO

    # set up the species annotation table
    spi = iface.sp_info
    spi.connect_after('size-allocate', lambda *args: [col.
                      set_max_width(spi.get_allocated_width() // 2) for col in spi.get_columns()])
    spi.set_model(data.sp_model)
    spi.set_reorderable(True)
    for i, title in enumerate(['sample', 'pid', 'species', 'runner-up']):
        crt = Gtk.CellRendererText(editable=True)
        crt.connect('edited', _save_sp_edit, spi, i, gui)
        col = Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=i)
        col.set_sort_column_id(i)
        col.set_resizable(True)
        spi.append_column(col)
    spi.columns_autosize()

    iface.blast_seen = False


def prep0(gui, path):
    data, iface = gui.data, gui.iface
    if not path:
        return False
    iface.i = 0
    iface.k = 3
    iface.text = 'search for BLAST+ installation'
    path = Path(path).resolve()
    if path.is_file():
        path = path.parent  # path is supposed to be a directory!
    data.blast_path = path  # save path to executable in project

    # check if all necessary BLAST+ scripts are available
    missing = [a for a in ['blastn', 'blastp', 'blastdbcmd', 'update_blastdb.pl',
                           'makeblastdb'] if not shutil.which(path / a)]
    if missing:
        shared.show_notification(gui, msg='Necessary scripts could not be '
                                          'found on the $PATH:', items=missing)
        sleep(.1)
        GObject.idle_add(prep4, gui)
        return True

    iface.blast_exe.set_filename(str(path))
    mo = Gtk.ListStore(str, str, str, str, bool)
    iface.blast_db.set_model(mo)  # shared model
    iface.db_info.set_model(mo)
    iface.db_info.set_reorderable(True)
    for i, title in enumerate(['database', 'size', 'title']):
        col = Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(),
                                 text=i + 1, underline=4)
        col.set_resizable(True)
        iface.db_info.append_column(col)

    iface.thread = threading.Thread(target=prep1, args=[gui, path, mo])
    GObject.timeout_add(50, shared.update, iface, PAGE)
    iface.i = 1
    iface.thread.start()
    return True


def prep1(gui, path, mo):
    """
    There is a BLAST+ executable on the $PATH. Pre-select it, search for local and remote
    databases with BLAST+ commands, then set-up the db_info and the blast_db ComboBox.
    :param gui:
    :param binary:
    :return:
    """
    data, iface = gui.data, gui.iface
    try:
        iface.text = 'look for local databases'
        output = subprocess.check_output(
            str(path / 'blastdbcmd') + ' -recursive -list_outfmt "%f\t%t\t%U" -list ' +
            str(BASE_DIR), shell=True).decode('utf-8').strip().split('\n')
        if output and output != ['']:
            sleep(.05)
            GObject.idle_add(prep2, iface.db_info, output, mo)
            iface.i = 2
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)

    try:
        iface.text = 'look for online pre-compiled databases'
        output = subprocess.check_output(str(path / 'update_blastdb.pl') +
                                         ' --source gcp --showall tsv',
                                         shell=True).decode('utf-8')
        sleep(.05)
        GObject.idle_add(prep3, data.remote_dbs, output, mo, iface.remote_db)
        iface.i = 3
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)
    sleep(.05)
    GObject.idle_add(prep4, gui)
    return True


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


def prep4(gui):
    iface = gui.iface
    iface.thread.join()
    # sync db selection in db_info TreeView and blast_db ComboBox
    sel = iface.db_info.get_selection()
    sel.connect('changed', lambda *args: iface.blast_db.
                set_active(sel.get_selected_rows()[1][0].get_indices()[0]))
    sel.select_path(0)
    iface.blast_db.connect('changed', lambda *args: sel
                           .select_path(iface.blast_db.get_active()))
    LOG.info('BLAST prep done')
    return False


def refresh(gui):
    """Re-view the page. Re-load the current genes to the marker-gene selector"""
    data, iface = gui.data, gui.iface

    if not iface.blast_seen:
        # look for BLAST+ executable on the $PATH and in the dataset; fill db_info table
        for path in [shutil.which('blastn'), data.blast_path]:
            if prep0(gui, path):
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
    shared.set_changed(gui, PAGE, False)


def start_BLAST(widget, gui, mode, *args):
    """Set-up the BLAST thread"""
    data, iface = gui.data, gui.iface

    if 'blast_wrapper' in iface:
        iface.pill2kill.set()
        if iface.blaster.is_alive():
            iface.blaster.stop()
        return

    iface.blast_spinner.start()
    iface.spinbox.set_visible(True)
    gene = iface.blast_gene.get_active_text()
    if mode == 'local':
        fasta = gui.wd / gene / (gene + '.fasta')
        if not fasta.is_file():
            shared.show_notification(gui, '%s with sequences for BLAST does not exist' % fasta)
            iface.blast_spinner.stop()
            iface.spinbox.set_visible(False)
            return
        pars = [True, False, False]
        adapt_button = iface.blast_local, mode
        # re-read all sequence data
        seqdata = {g: {r.id: r for r in SeqIO.parse(
            gui.wd / g / (g + '_all.fasta'), 'fasta')} for g in data.genes}

    elif mode == 'remote':
        fasta = iface.missing_fasta.get_filename()
        if not fasta and mode == 'remote':
            shared.show_notification(gui, 'select a FASTA with sequence data first')
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
            gui.wd / g / (g + '_all.fasta'), 'fasta')} for g in data.genes}

    pars = dict(zip(['no_remote', 'no_local', 'BLAST_xml'], pars))
    (im, la), tx = adapt_button[0].get_child().get_children(), adapt_button[1]
    im.set_from_icon_name('media-playback-stop-symbolic', 4)
    la.set_text('Stop')
    db_row = iface.blast_db.get_model()[iface.blast_db.get_active()]

    # define static parameters
    args = {'gui': True,
            'no_BLAST': False,
            'genes': [gene],
            'db': db_row[1],
            'remote_db': data.remote_dbs[iface.remote_db.get_active()][0],
            'timeout': DOWNLOAD_TIMEOUT,
            'dir': gui.wd,
            'xml': gui.wd / PATHS.xml,
            'www_xml': str(gui.wd / PATHS.www_xml),
            'tsv': gui.wd / PATHS.tsv,
            'bad_seqs': gui.wd / PATHS.bad_seqs,
            'missing_fasta': gui.wd / PATHS.missing_fasta,
            'dbpath': Path(db_row[0]).parent if db_row[0] else None}
    args.update(pars)
    args = Namespace(**args)
    reader = Namespace(**{'seqdata': seqdata, 'metadata': None})

    iface.tup = (im, la, tx, gene)
    iface.blaster = blast.blast_build(args, reader)
    iface.blast_wrapper = threading.Thread(
        # iface.blaster = multiprocessing.Process(
        target=do_BLAST, args=[gui])
    iface.blast_wrapper.start()
    return
    # return to main loop


def do_BLAST(gui):
    """Run BLAST thread"""
    with shared.bump_log_level(LOG):
        gui.iface.blaster.start()
        gui.iface.blaster.join()
    GObject.idle_add(stop_BLAST, gui)


def stop_BLAST(gui):
    """Finish the BLAST thread"""
    data, iface = gui.data, gui.iface
    if 'blaster' in iface:
        if iface.blaster.update:  # a local database was created
            iface.blast_db.get_model()[iface.blast_db.get_active()][-1] = True
        if iface.blaster.missing_fasta.is_file():
            iface.missing_fasta.set_filename(str(iface.blaster.missing_fasta))

    im, la, tx, gene = iface.tup

    _refill(gui, gene, fresh=True)

    iface.blast_spinner.stop()
    iface.spinbox.set_visible(False)
    im.set_from_icon_name('media-playback-start-symbolic', 4)
    la.set_text(tx)

    if iface.pill2kill.is_set():
        tx = 'stopped BLAST'
        shared.show_notification(gui, msg=tx, stay_secs=2)
    else:
        tx = 'BLAST finished'
        iface.blast_help_stack.set_visible_child_name('sp_info')
        if iface.notebook.get_current_page() != PAGE:
            shared.show_notification(gui, msg=tx, stay_secs=5)
    LOG.info(tx)
    iface.pill2kill.clear()
    del iface.blast_wrapper, iface.blaster
    shared.set_changed(gui, PAGE, False)


def _save_custom_remote_db(entry, *args):
    combo = entry.get_parent().get_parent()
    if not combo.get_active_iter():
        tx = entry.get_text()
        combo.get_model().append([tx, -1])
        LOG.debug('entered remote_db %s' % tx)


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


def _refill(gui, gene=None, fresh=False):
    """Re-fill the species annotation table"""
    data, iface = gui.data, gui.iface
    LOG.debug('re-fill species annotations')
    if not gene:
        gene = iface.blast_gene.get_active_text()
    if 'df' not in iface.tempspace or fresh:
        iface.tempspace.df = _df_with_sp(gui.wd / PATHS.tsv)
    df = iface.tempspace.df
    data.sp_model.clear()
    df = df.loc[df['gene'] == gene]
    for sample, r in df.iterrows():
        if 'pid' not in r or isnan(r.pid):  # TODO errors ?
            pid, r.BLAST_species, r.extra_species = '', '', ''
        elif type(r.pid) == float:
            pid = '%.2f' % r.pid if r.pid < 100 else '100'
        else:
            pid = r.pid
        data.sp_model.append([sample, pid, r.BLAST_species, r.extra_species])


def _save_sp_edit(cell, path, new_text, tv, col, gui):
    data, iface = gui.data, gui.iface
    if 'df' not in gui.iface.tempspace:
        gui.iface.tempspace.df = _df_with_sp(gui.wd / PATHS.tsv)
    df = gui.iface.tempspace.df
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
    shared.save_row_edits(cell, path, str(new_text), tv, col)
    shared.set_changed(gui, PAGE)


def human_bytes(num):
    step_unit = 1024.0
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < step_unit:
            return '{:.0f} {}'.format(num, x)
        num /= step_unit
