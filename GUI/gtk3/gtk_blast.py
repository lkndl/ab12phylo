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

    # iface.blast_spinner.props.visible = False
    # iface.blast_spinner.props.active = False
    # iface.blast_spinner.stop()

    # connect buttons
    iface.blast_local.connect('clicked', start_BLAST, gui, 'local')
    iface.blast_remote.connect('clicked', start_BLAST, gui, 'remote')
    iface.xml_import.connect('file_set', start_BLAST, gui, [iface.xml_import.get_filename()])

    # look for BLAST+ executable on the $PATH
    binary = shutil.which('blastn')
    if binary:
        iface.thread = threading.Thread(target=do_prep, args=[gui, binary])
        iface.running = True
        iface.thread.start()

    # TODO iface.remote_db  # ComboBoxText: editable. handle db error?


def do_prep(gui, binary):
    """
    There is a BLAST+ executable on the $PATH. Pre-select it, search for local and remote
    databases with BLAST+ commands, then set-up the db_info and the blast_db ComboBox.
    :param gui:
    :param binary:
    :return:
    """
    data, iface = gui.data, gui.iface

    iface.blast_exe.set_filename(binary)
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
    # iface.blast_db.set_entry_text_column(0)
    crt = Gtk.CellRendererText()
    iface.blast_db.pack_start(crt, True)
    iface.blast_db.add_attribute(crt, 'text', 1)
    iface.blast_db.add_attribute(crt, 'underline', 4)
    tv.columns_autosize()

    # sync db selection in db_info TreeView and blast_db ComboBox
    sel = iface.db_info.get_selection()
    iface.blast_db.connect('changed', lambda args: sel.select_path(iface.blast_db.get_active()))
    sel.connect('changed', lambda args: iface.blast_db.set_active(sel.get_selected_rows()[1][0].get_indices()[0]))

    # iface.blast_db.set_entry_text_column(0)
    # [iface.blast_db.append_text(db) for db in shared.get_column(mo, 1)]
    # iface.blast_db.set_active(0)
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
    elif not shared.get_changed(gui, PAGE):
        shared.show_notification(gui, 'BLAST already finished, please proceed')
        return
    # TODO early aborts?
    # elif all([(gui.wd / gene / ('%s_raw_msa.fasta' % gene)).exists() for gene in data.genes]) \
    #         and run_after and not shared.get_errors(gui, PAGE):  # b)
    #     shared.set_changed(gui, PAGE, False)
    #     [do_func(gui) for do_func in run_after]
    #     return

    # iface.blast_spinner.props.visible = True
    # iface.blast_spinner.props.active = True
    # iface.blast_spinner.start()

    # define mode-dependant parameters
    pars = [True, True, mode]
    if type(mode) == str:
        if mode == 'local':
            pars = [False, False, False]
        elif mode == 'remote':
            pars = [False, True, False]
    pars = dict(zip(['no_remote', 'no_local', 'BLAST_xml'], pars))

    db_row = iface.blast_db.get_model()[iface.blast_db.get_active()]

    # define static parameters
    args = {'no_BLAST': False,
            'genes': [gene],
            'db': db_row[1],
            'remote_db': iface.remote_db.get_active_text(),
            'timeout': 20,
            'dir': gui.wd,
            'xml': gui.wd / PATHS.xml,
            'www_xml': str(gui.wd / PATHS.www_xml),
            'tsv': gui.wd / PATHS.tsv,
            'bad_seqs': None,
            'missing_fasta': gui.wd / PATHS.missing_fasta,
            'dbpath': Path(db_row[0]).parent}
    args.update(pars)
    args = Namespace(**args)

    records = {record.id: record for record in SeqIO.parse(fasta, 'fasta')}

    reader = Namespace(**{'seqdata': {gene: records},
                          'metadata': {gene: {k: {} for k in records.keys()}}})
    # TODO repair with actual metadata.tsv

    # MARK thread in a thread seems to work?
    iface.blaster = threading.Thread(target=do_BLAST, args=[gui, blast.blast_build(args, reader)])
    iface.blast_spinner.start()
    iface.blaster.start()
    return
    # return to main loop



def do_BLAST(gui, blaster):
    """Run BLAST thread."""
    blaster.start()

    GObject.idle_add(stop_BLAST, gui)


def stop_BLAST(gui):
    """Finish the Gblocks thread"""
    data, iface = gui.data, gui.iface
    iface.blaster.join()
    # iface.blast_spinner.props.visible = False
    # iface.blast_spinner.props.active = False
    iface.blast_spinner.stop()
    gui.win.show_all()
    LOG.info('BLAST finished')
    shared.set_changed(gui, PAGE, False)


def human_bytes(num):
    step_unit = 1024.0
    for x in ['B', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < step_unit:
            return '{:.0f} {}'.format(num, x)
        num /= step_unit
