# 2020 Leo Kaindl

import logging

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
        self.genes = list()  # used also before seqdata exists
        self.csvs = dict()
        self.seqdata = dict()
        self.metadata = dict()
        self.seed = 0
        self.record_order = list()
        self.qal_model = picklable_liststore(str,  # id
                                             int,  # has phreds
                                             bool)  # low quality
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
        self.gene_ids = dict()
        self.gene_for_preview = ''
        self.sp_model = picklable_liststore(str,  # id
                                            str,  # pid
                                            str,  # species
                                            str)  # extra_species
        self.blast_path = None  # for non-$PATH BLAST+ executable

    def new_project(self):
        self.overwrite(project_dataset())

    def overwrite(self, new_dataset):
        for attr in [a for a in dir(self) if not callable(getattr(self, a)) and not a.startswith('__')]:
            old = self.__getattribute__(attr)
            if type(old) == picklable_liststore:
                old.clear()
                try:
                    [old.append(row[:]) for row in new_dataset.__getattribute__(attr)]
                except AttributeError:
                    pass  # try to go without this model
            elif type(old) == dict:
                old.clear()
                try:
                    old.update(new_dataset.__getattribute__(attr))
                except AttributeError:
                    pass  # try to go without dict values
            else:
                try:
                    self.__setattr__(attr, new_dataset.__getattribute__(attr))
                except (AttributeError, TypeError) as ex:
                    LOG.error(ex)


class picklable_liststore(Gtk.ListStore):
    """
    kudos go to samplebias on https://stackoverflow.com/a/5969700
    """

    def __reduce__(self):
        rows = list()
        try:
            rows = [list(row) for row in self]
            coltypes = [type(c) for c in rows[0]]
            return _unpickle_liststore, (self.__class__, coltypes, rows)
        except IndexError as ex:
            cols = self.get_n_columns()
            # allow saving of empty data stores
            if cols == 1:
                coltypes = [str]
            if cols == 2:
                coltypes = [str, str]
            elif cols == 4:
                coltypes = [str, str, str, str]
            elif cols == 3:
                coltypes = [str, int, bool]
            elif cols == 7:
                coltypes = [str, str, str, str, str, bool, bool, str]
            return _unpickle_liststore, (self.__class__, coltypes, rows)
        except Exception as ex:
            LOG.exception(ex)


def _unpickle_liststore(cls, col_types, rows):
    inst = cls.__new__(cls)
    inst.__init__(*col_types)
    for row in rows:
        inst.append(row)
    return inst