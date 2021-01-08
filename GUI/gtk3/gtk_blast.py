# 2020 Leo Kaindl

import logging
import shutil
import subprocess
import threading
from argparse import Namespace
from pathlib import Path

import gi
import pandas as pd
from Bio import SeqIO

from ab12phylo import blast

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from GUI.gtk3 import shared
from static import PATHS, BASE_DIR

LOG = logging.getLogger(__name__)
PAGE = 5


def init(gui):
    """Initialize the page. Connect buttons."""
    data, iface = gui.data, gui.iface

    iface.blast_spinner.stop()
    iface.spinbox.set_visible(False)

    # connect buttons
    iface.blast_local.connect('clicked', start_BLAST, gui, 'local')
    iface.blast_remote.connect('clicked', start_BLAST, gui, 'remote')
    iface.xml_import.connect('file_set', start_BLAST, gui, [iface.xml_import.get_filename()])

    iface.blast_exe.connect('file_set', lambda args: start_prep(gui, iface.blast_exe.get_filename()))

    # look for BLAST+ executable on the $PATH
    start_prep(gui, shutil.which('blastn'))

    # set up the species annotation table
    iface.sp_info.set_model(data.sp_model)
    for i, title in enumerate(['sample', 'pid', 'hits', 'species', 'others']):
        iface.sp_info.append_column(
            Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(), text=i))
        iface.sp_info.get_column(i).set_sort_column_id(i)
    iface.sp_info.columns_autosize()

    # TODO iface.remote_db  # ComboBoxText: editable. handle db error?
    # TODO self.blast_path = None  # for non-$PATH BLAST+ executable
    # TODO Test for import-only!


def start_prep(gui, path):
    if not path:
        return
    path = Path(path).resolve()
    if path.is_file():
        path = path.parent
    gui.iface.thread = threading.Thread(target=do_prep, args=[gui, path])
    gui.iface.running = True
    gui.iface.thread.start()


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
        GObject.idle_add(stop_prep, gui)
        return

    iface.blast_exe.set_filename(str(path))
    mo = Gtk.ListStore(str, str, str, str, bool)
    tv = iface.db_info
    tv.set_model(mo)
    for i, title in enumerate(['database', 'size', 'title']):
        tv.append_column(
            Gtk.TreeViewColumn(title=title, cell_renderer=Gtk.CellRendererText(), text=i + 1, underline=4))

    # look for locally installed databases
    try:
        output = subprocess.check_output('blastdbcmd -recursive -list_outfmt "%f\t%t\t%U" -list ' + str(BASE_DIR),
                                         shell=True).decode('utf-8').strip().split('\n')
        if output and output != ['']:
            tv.set_tooltip_text('Underlining indicates a local database')
            for db in output:
                db = db.split('\t')
                try:
                    mo.append([db[0], Path(db[0]).name, human_bytes(
                        sum(f.stat().st_size for f in Path(db[0]).parent.glob('**/*') if f.is_file())), db[1], True])
                except IndexError:
                    pass  # if directory does not exist or is empty
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)

    # look for only pre-compiled databases
    try:
        output = subprocess.check_output('update_blastdb.pl --source gcp --showall tsv',
                                         shell=True).decode('utf-8')
        for db in output.strip().split('\n')[1:]:
            db = db.split('\t')
            mo.append([None, db[0], human_bytes(2 ** 30 * float(db[2].strip())), db[1], False])
    except subprocess.CalledProcessError as ex:
        LOG.error(ex)

    # setup blast_db selection ComboBox.
    iface.blast_db.set_model(mo)
    crt = Gtk.CellRendererText()
    iface.blast_db.pack_start(crt, True)
    iface.blast_db.add_attribute(crt, 'text', 1)
    iface.blast_db.add_attribute(crt, 'underline', 4)
    # iface.blast_db.set_active(0)
    tv.columns_autosize()

    # sync db selection in db_info TreeView and blast_db ComboBox
    sel = iface.db_info.get_selection()
    sel.connect('changed', lambda args: iface.blast_db.set_active(sel.get_selected_rows()[1][0].get_indices()[0]))
    sel.select_path(0)
    iface.blast_db.connect('changed', lambda args: sel.select_path(iface.blast_db.get_active()))

    GObject.idle_add(stop_prep, gui)


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

    # allow skipping BLAST
    shared.set_changed(gui, PAGE, False)


def start_BLAST(widget, gui, mode, *args):
    """Set-up the BLAST thread."""
    data, iface = gui.data, gui.iface

    gene = gui.iface.blast_gene.get_active_text()
    fasta = gui.wd / gene / (gene + '.fasta')
    if not fasta.is_file():
        shared.show_notification(gui, '%s with sequences for BLAST does not exist' % fasta)
        return

    if 'blaster' in iface:
        shared.show_notification(gui, 'BLAST already running')
        return
    # elif not shared.get_changed(gui, PAGE):
    #     shared.show_notification(gui, 'BLAST already finished, please proceed')
    #     return
    # TODO early aborts?
    # elif all([(gui.wd / gene / ('%s_raw_msa.fasta' % gene)).exists() for gene in data.genes]) \
    #         and run_after and not shared.get_errors(gui, PAGE):  # b)
    #     shared.set_changed(gui, PAGE, False)
    #     [do_func(gui) for do_func in run_after]
    #     return

    iface.blast_spinner.start()
    iface.spinbox.set_visible(True)

    # convert metadata dict do pandas DataFrame
    df = pd.concat({gene: pd.DataFrame.from_dict(
        data.metadata[gene], orient='index') for gene in data.metadata.keys()})
    df.index.names = ['gene', 'id']
    df.reset_index(level=0, inplace=True)
    df.to_csv(gui.wd / PATHS.tsv, sep='\t', na_rep='', header=True, index=True)

    # define mode-dependant parameters
    pars = [True, True, mode]
    adapt_button = None
    if type(mode) == str:
        if mode == 'local':
            pars = [True, False, False]
            adapt_button = iface.blast_local
        elif mode == 'remote':
            pars = [False, True, False]
            adapt_button = iface.blast_remote
    pars = dict(zip(['no_remote', 'no_local', 'BLAST_xml'], pars))

    if adapt_button:
        im, la = adapt_button.get_child().get_children()
        im.set_from_icon_name('media-playback-stop-symbolic', 4)
        la.set_text('Stop')

    db_row = iface.blast_db.get_model()[iface.blast_db.get_active()]

    # define static parameters
    args = {'gui': True,
            'no_BLAST': False,
            'genes': [gene],
            'db': db_row[1],
            'remote_db': iface.remote_db.get_active_text(),
            'timeout': 20,
            'dir': gui.wd,
            'xml': gui.wd / PATHS.xml,
            'www_xml': str(gui.wd / PATHS.www_xml),
            'tsv': gui.wd / PATHS.tsv,
            'bad_seqs': gui.wd / PATHS.bad_seqs,
            'missing_fasta': gui.wd / PATHS.missing_fasta,
            'dbpath': Path(db_row[0]).parent}
    args.update(pars)
    args = Namespace(**args)

    records = {record.id: record for record in SeqIO.parse(fasta, 'fasta')}
    reader = Namespace(**{'seqdata': {gene: records}, 'metadata': None})
    iface.blaster = threading.Thread(
        target=do_BLAST, args=[gui, blast.blast_build(args, reader),
                               (im, la, mode) if adapt_button else None])
    iface.blaster.start()
    return
    # return to main loop


def do_BLAST(gui, blaster, tup):
    """Run BLAST thread."""
    blaster.start()
    blaster.join()
    # fill the table
    blaster.df.reset_index('gene', inplace=True)
    for sample, r in blaster.df.loc[blaster.df['gene'] == blaster.gene].iterrows():
        pid = '' if pd.isna(r.pid) else '%.2f' % r.pid if r.pid < 100 else '100'
        gui.data.sp_model.append([sample, pid, r.hit_ratio, r.BLAST_species, r.extra_species])
    # flip the stack
    gui.iface.blast_help_stack.set_visible_child_name('sp_info')
    GObject.idle_add(stop_BLAST, gui, tup)


def stop_BLAST(gui, tup):
    """Finish the Gblocks thread"""
    data, iface = gui.data, gui.iface
    iface.blast_spinner.stop()
    iface.spinbox.set_visible(False)
    if tup:
        tup[0].set_from_icon_name('media-playback-start-symbolic', 4)
        tup[1].set_text('Run' if tup[2] == 'local' else 'Submit')
    del iface.blaster
    LOG.info('BLAST finished')
    shared.set_changed(gui, PAGE, False)


def human_bytes(num):
    step_unit = 1024.0
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < step_unit:
            return '{:.0f} {}'.format(num, x)
        num /= step_unit


def pd_from_(metadata):
    # convert metadata dict do pandas DataFrame
    df = pd.concat({gene: pd.DataFrame.from_dict(
        metadata[gene], orient='index') for gene in metadata.keys()})
    df.index.names = ['gene', 'id']
    df.reset_index(level=0, inplace=True)
    return df
