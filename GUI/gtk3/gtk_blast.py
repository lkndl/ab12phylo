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
from numpy import isnan

from ab12phylo import blast

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import shared
from static import PATHS, BASE_DIR, DOWNLOAD_TIMEOUT

LOG = logging.getLogger(__name__)
PAGE = 5


def init(gui):
    """Initialize the page. Connect buttons."""
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
    iface.blast_exe.connect('file_set', lambda *args: start_prep(gui, iface.blast_exe.get_filename()))
    iface.blast_gene.connect('button-release-event', lambda co, *args: co.
                             popdown() if co.props.popup_shown else co.popup())

    # set up the species annotation table
    spi = iface.sp_info
    spi.connect_after('size-allocate', lambda *args: [col.
                      set_max_width(spi.get_allocated_width() // 2) for col in spi.get_columns()])
    spi.set_model(data.sp_model)
    spi.set_reorderable(True)
    for i, title in enumerate(['sample', 'pid', 'species', 'runner-up']):
        crt = Gtk.CellRendererText(editable=True)
        crt.connect('edited', shared.save_row_edits, spi, i)
        col = Gtk.TreeViewColumn(title=title, cell_renderer=crt, text=i)
        col.set_sort_column_id(i)
        col.set_resizable(True)
        spi.append_column(col)
    spi.columns_autosize()

    # look for BLAST+ executable on the $PATH and in the dataset; fill db_info table
    for path in [shutil.which('blastn'), data.blast_path]:
        if start_prep(gui, path):
            break


def start_prep(gui, path):
    data, iface = gui.data, gui.iface
    if not path:
        return False
    path = Path(path).resolve()
    if path.is_file():
        path = path.parent  # path is supposed to be a directory!
    gui.data.blast_path = path
    iface.thread = threading.Thread(target=do_prep, args=[gui, path])
    iface.running = True
    iface.thread.start()
    return True


def do_prep(gui, path):
    """
    There is a BLAST+ executable on the $PATH. Pre-select it, search for local and remote
    databases with BLAST+ commands, then set-up the db_info and the blast_db ComboBox.
    :param gui:
    :param binary:
    :return:
    """
    data, iface = gui.data, gui.iface

    # check if all necessary BLAST+ scripts are available
    missing = [a for a in ['blastn', 'blastp', 'blastdbcmd', 'update_blastdb.pl',
                           'makeblastdb'] if not shutil.which(path / a)]
    if missing:
        shared.show_notification(gui, msg='Necessary scripts could not be '
                                          'found on the $PATH:', items=missing)
        sleep(.1)
        GObject.idle_add(stop_prep, gui)
        return True

    iface.blast_exe.set_filename(str(path))
    mo = Gtk.ListStore(str, str, str, str, bool)
    tv = iface.db_info
    tv.set_model(mo)
    tv.set_reorderable(True)
    for i, title in enumerate(['database', 'size', 'title']):
        col = Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(),
                                 text=i + 1, underline=4)
        col.set_resizable(True)
        tv.append_column(col)

    # look for locally installed databases
    try:
        output = subprocess.check_output(
            'blastdbcmd -recursive -list_outfmt "%f\t%t\t%U" -list ' + str(BASE_DIR),
            shell=True).decode('utf-8').strip().split('\n')
        if output and output != ['']:
            tv.set_tooltip_text('Underlining indicates a local database')
            for db in output:
                db = db.split('\t')
                try:
                    mo.append([db[0], Path(db[0]).name, human_bytes(
                        sum(f.stat().st_size for f in Path(db[0]).parent
                            .glob('**/*') if f.is_file())), db[1], True])
                except IndexError:
                    pass  # if directory does not exist or is empty
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)

    # look for online pre-compiled databases
    try:
        output = subprocess.check_output('update_blastdb.pl --source gcp --showall tsv',
                                         shell=True).decode('utf-8')
        data.remote_dbs.clear()
        for i, db in enumerate(output.strip().split('\n')[1:]):
            db = db.split('\t')
            mo.append([None, db[0], human_bytes(
                2 ** 30 * float(db[2].strip())), db[1], False])
            # also fill remote_db
            data.remote_dbs.append([db[0], i])
            if db[0] == 'nt':
                iface.remote_db.set_active(i)
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)

    # setup blast_db selection ComboBox
    iface.blast_db.set_model(mo)

    # sync db selection in db_info TreeView and blast_db ComboBox
    sel = iface.db_info.get_selection()
    sel.connect('changed', lambda *args: iface.blast_db.set_active(
        sel.get_selected_rows()[1][0].get_indices()[0]))
    sel.select_path(0)
    iface.blast_db.connect('changed', lambda *args: sel.select_path(iface.blast_db.get_active()))

    sleep(.1)
    GObject.idle_add(stop_prep, gui)
    return True


def stop_prep(gui):
    iface = gui.iface
    iface.running = False
    iface.thread.join()
    LOG.info('BLAST prep done')
    return


def refresh(gui):
    """Re-view the page. Re-load the current genes to the marker-gene selector"""
    data, iface = gui.data, gui.iface

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
    """ Set-up the BLAST thread. """
    data, iface = gui.data, gui.iface

    if 'blast_wrapper' in iface:
        iface.pill2kill.set()
        if iface.blaster.is_alive():
            iface.blaster.stop()
        shared.show_notification(gui, 'stopped BLAST%s' %
                                 (':\nremote might not die' if mode == 'remote' else ''))
        return

    iface.blast_spinner.start()
    iface.spinbox.set_visible(True)
    gene = gui.iface.blast_gene.get_active_text()
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

    iface.blaster = blast.blast_build(args, reader)
    iface.blast_wrapper = threading.Thread(
        # iface.blaster = multiprocessing.Process(
        target=do_BLAST, args=[gui, (im, la, tx, gene)])
    iface.blast_wrapper.start()
    return
    # return to main loop


def do_BLAST(gui, tup):
    """Run BLAST thread"""
    blaster = gui.iface.blaster
    pill2kill = gui.iface.pill2kill
    with shared.bump_log_level(LOG):
        while not pill2kill.wait(.5):
            blaster.start()
            blaster.join()
            blaster.stop()
            pill2kill.set()

    blaster.stop()
    GObject.idle_add(stop_BLAST, gui, tup)


def stop_BLAST(gui, tup):
    """Finish the BLAST thread"""
    data, iface = gui.data, gui.iface
    if gui.iface.blaster.update:  # a local database was created
        iface.blast_db.get_model()[iface.blast_db.get_active()][-1] = True
    del iface.blast_wrapper, gui.iface.blaster
    iface.pill2kill.clear()
    im, la, tx, gene = tup

    # re-fill the table
    gui.data.sp_model.clear()
    df = pd.read_csv(gui.wd / PATHS.tsv, sep='\t', dtype={'id': str})
    df = df.set_index('id')
    df = df.loc[df['gene'] == gene]
    df.extra_species.fillna('', inplace=True)
    df.BLAST_species.fillna('', inplace=True)
    for sample, r in df.iterrows():
        pid = '' if isnan(r.pid) else '%.2f' % r.pid if r.pid < 100 else '100'
        gui.data.sp_model.append([sample, pid, r.BLAST_species, r.extra_species])
    gui.iface.blast_help_stack.set_visible_child_name('sp_info')

    iface.blast_spinner.stop()
    iface.spinbox.set_visible(False)
    im.set_from_icon_name('media-playback-start-symbolic', 4)
    la.set_text(tx)
    LOG.info('BLAST finished')
    shared.set_changed(gui, PAGE, False)
    if iface.notebook.get_current_page() != PAGE:
        shared.show_notification(gui, msg='BLAST finished')


def _save_custom_remote_db(entry, *args):
    combo = entry.get_parent().get_parent()
    if not combo.get_active_iter():
        tx = entry.get_text()
        combo.get_model().append([tx, -1, True])
        LOG.debug('entered remote_db %s' % tx)


def human_bytes(num):
    step_unit = 1024.0
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < step_unit:
            return '{:.0f} {}'.format(num, x)
        num /= step_unit


def pd_from_(metadata):  # TODO delete
    # convert metadata dict do pandas DataFrame
    df = pd.concat({gene: pd.DataFrame.from_dict(
        metadata[gene], orient='index') for gene in metadata.keys()})
    df.index.names = ['gene', 'id']
    df.reset_index(level=0, inplace=True)
    return df
