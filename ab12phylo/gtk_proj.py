# 2021 Leo Kaindl
import argparse
import logging
from argparse import Namespace

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

LOG = logging.getLogger(__name__)


class project_dataset:
    def __init__(self):
        self.trace_store = picklable_liststore(str,  # path
                                               str,  # filename
                                               str,  # well/id
                                               str,  # plate
                                               str,  # gene
                                               bool,  # reference
                                               bool,  # reversed
                                               str)  # color

        self.plate_store = picklable_liststore(str,  # path
                                               str,  # filename
                                               str,  # plate ID
                                               str)  # errors
        self.rx_fired = [False] * 12
        # rx_fired stores if these columns have been parsed before:
        # [bool] with trace and then plate columns. sum(iface.rx_fired) = 5
        self.genes = list()  # used also before seqdata exists
        self.csvs = dict()
        self.seqdata = dict()
        self.metadata = dict()
        self.seed = 0
        self.record_order = list()
        self.qal_model = picklable_liststore(str,  # id
                                             str,  # gene
                                             int,  # has phreds
                                             bool)  # low quality
        self.search_rev = False
        self.rgx = Namespace()
        self.qal = Namespace(gene_roll='all', accept_rev=True, accept_nophred=True,
                             min_phred=30, trim_out=8, trim_of=10, bad_stretch=5)
        self.msa = Namespace(cmd=dict(), remote_cmd=dict(), last_cmd='')
        self.gbl = Namespace(ignore_ids=set())
        self.ml = Namespace(evo_model='GTR', evo_modify='+G4', cpu_count='1', cpu_use='1',
                            tool='raxml-ng', ml_cmd='command not yet defined')
        self.phy = Namespace(gap_share=.1, unk_share=.1, flip=.7, dist=0, space=1, sel_gene=None,
                             did_BLAST=False, query='', exclude='',
                             rect=True, circ=False, unro=False, tbe=True, fbp=False, supp=True,
                             spec=True, axis=True, align=False, pmsa=False,
                             pdf=False, svg=False, png=True, nwk=True, html=False,
                             h=dict(), qh=dict(), modify=dict())

        self.gbl_model = picklable_liststore(str)  # id
        # set up indicator of changes, tabs are not disabled initially
        self.change_indicator = [False] * 20
        self.errors_indicator = [False] * 20
        self.page = 0
        self.qal_shape = [0, 0]
        self.gbl_shape = [0, 0]  # width-height
        self.msa_shape = [0, 0, 0, 0]  # width-height before and after trimming
        self.msa_lens = list()
        self.msa_hash = ''
        self.gene_ids = dict()  # {gene: {seq_ids}}
        self.gene_for_preview = ''
        self.sp_model = picklable_liststore(str,  # id
                                            str,  # pid
                                            str,  # species
                                            str,  # extra_species
                                            str)  # low quality indicator color
        self.remote_dbs = picklable_liststore(str,  # name
                                              int)  # for synchronization with ComboBox, i think ...
        self.blast_path = None  # for non-$PATH BLAST+ executable
        self.evo_models = picklable_liststore(str, str)  # id, use_file
        [self.evo_models.append([evo, None]) for evo in ['GTR', 'JC', 'SYM', 'K80', 'K81', 'HKY', 'TN93',
                                                         'TPM2', 'TPM3', 'TIM2', 'TIM3', 'TVM', 'from file']]
        self.pop_model = picklable_liststore(str,  # gene
                                             int,  # valid sites
                                             int,  # S
                                             float,  # k
                                             float,  # pi
                                             float,  # theta_w
                                             float,  # D_T
                                             int,  # unique seqs
                                             int,  # gaps
                                             int)  # unknown
        self.tree_anno_model = picklable_liststore(str,  # sample ID
                                                   str,  # species
                                                   str,  # foreground color
                                                   str)  # background color
        self.colors = None

    def new_project(self):
        self.overwrite(project_dataset())

    def overwrite(self, new_dataset):
        for attr in [a for a in dir(self) if not callable(getattr(self, a)) and not a.startswith('__')]:
            old = self.__getattribute__(attr)
            if type(old) == picklable_liststore:
                old.clear()
                try:
                    [old.append(row[:]) for row in new_dataset.__getattribute__(attr)]
                except (AttributeError, ValueError):
                    pass  # try to go without this model
            elif type(old) == dict:
                old.clear()
                try:
                    old.update(new_dataset.__getattribute__(attr))
                except (TypeError, AttributeError):
                    pass  # try to go without dict values
            elif type(old) == argparse.Namespace:
                try:
                    for k, v in new_dataset.__getattribute__(attr).__dict__.items():
                        old.__setattr__(k, v)
                except (AttributeError, TypeError, ValueError):
                    pass  # let's be courageous
            else:
                try:
                    self.__setattr__(attr, new_dataset.__getattribute__(attr))
                except (AttributeError, TypeError) as ex:
                    LOG.error(ex)


class picklable_liststore(Gtk.ListStore):
    """kudos go to samplebias on https://stackoverflow.com/a/5969700"""

    def __init__(self, *args):
        Gtk.ListStore.__init__(self, *args)
        self.col_types = args

    def __reduce__(self):
        rows = [list(row) for row in self]
        return _unpickle_liststore, (self.__class__, self.col_types, rows)

    def get_column(self, col_idx):
        """Extract a column from a Gtk.ListStore, because they are annoying"""
        col = list()
        for row in self:
            col.append(row[col_idx])
        return col


def _unpickle_liststore(cls, col_types, rows):
    inst = cls.__new__(cls)
    inst.__init__(*col_types)
    for row in rows:
        inst.append(row)
    return inst
