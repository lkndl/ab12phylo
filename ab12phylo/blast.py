# 2020 Leo Kaindl

"""
Runs a BLAST+ search on a local database and a remote NCBI BLAST search in the nucleotides
database for seqs missing from it.
Result :code:`XML` files are parsed, and dictionary of non-sequence data is updated and
written to file. Processes are run in a separate thread due to the Python GIL.
"""

import logging
import os
import shutil
import subprocess
import threading
from os import path
from time import time

import pandas
from numpy import nan
from Bio import SearchIO


class blast_build(threading.Thread):
    """
    Looks for pre-installed BLAST+. If 2.9.0 or newer, updates or downloads the database
    (usually ITS_RefSeq_Fungi). If no or an outdated BLAST+ installation is found,
    module exits with error.
    Updates and searches in local BLAST+ db with blastn and parses resulting .xml to the
    metadata .tsv. Runs an online NCBI BLAST for missing seqs.
    """

    def __init__(self, _args, reader):
        threading.Thread.__init__(self)
        self.log = logging.getLogger(__name__)

        # modes / variables
        self.no_BLAST = _args.no_BLAST
        if self.no_BLAST is True:
            return
        self.no_remote = _args.no_remote
        self.gene = _args.genes[0]
        self.db = _args.db
        self.remote_db = _args.remote_db
        self.timeout = _args.timeout

        # Files / Paths
        self.FASTA = path.join(_args.dir, self.gene, self.gene + '.fasta')
        self.XML = _args.xml
        self.www_XML = _args.www_xml
        self.TSV = _args.tsv
        self.bad_seqs = _args.bad_seqs
        self.missing_fasta = _args.missing_fasta
        self.dbpath = _args.dbpath
        self.blast_dir = None

        # data in memory
        self.seqdata = reader.seqdata
        self.metadata = reader.metadata

        # if BLAST db path was provided, use and don't update
        if self.dbpath is not None:
            self.update = False
        else:
            # save BLAST db in root directory
            self.dbpath = path.join(path.dirname(path.dirname(__file__)), 'blastdb', self.db)
            self.update = True
            os.makedirs(self.dbpath, exist_ok=True)

        try:
            binary = shutil.which('blastn')
            if binary is None:
                raise ValueError('BLAST+ not installed (not on the $PATH)')
            output = subprocess.check_output(binary + ' -version', shell=True).decode('utf-8')
            version = [int(i) for i in output.split(',')[0].split(' ')[-1].split('.')]
            if version[0] >= 2 and version[1] >= 9:
                self.blast_dir = path.dirname(binary)
            else:
                raise ValueError('BLAST+ in %s is outdated' % self.blast_dir)
        except (subprocess.CalledProcessError, ValueError) as e:
            self.log.error(e)
            exit(e)

    def run(self):
        # skip all BLASTing
        if self.no_BLAST is True:
            self.log.info('--SKIPPING BLAST--')  # leave it alone!
            return

        # run local BLAST
        not_found = self._blast_local()

        # if local BLAST failed, continue without annotation
        if not_found is None:
            return

        # online BLAST for seqs missing from local db
        if len(not_found) == 0:
            self.log.info('all seqs found in local %s' % self.db)
            return

        # save missing sequences to file
        with open(self.missing_fasta, 'w') as fh:
            fh.write('\n'.join([self.seqdata[self.gene][seq_id].format('fasta') for seq_id in not_found]))
            self.log.debug('sequences not found in %s saved in %s' % (self.db, self.missing_fasta))

        if self.no_remote is True:
            self.log.warning('skip online BLAST for missing seqs')
            # try reading old result
            if os.path.isfile(self.www_XML):
                try:
                    self.log.warning('reading XML result of previous online BLAST')
                    self._parse_remote_result()
                except Exception as ex:
                    self.log.error(ex)
            return

        self._blast_remote(not_found)

    def _blast_local(self):
        """
        Optionally update a local BLAST+ database, then check integrity.
        Run a BLAST+ search and finally write updated non-seq data to a .csv file.
        :return: missing sequences as list of sample IDs
        """
        # update database
        if self.update:
            self.log.debug('updating BLAST+ db ...')
            arg = path.join(self.blast_dir, 'update_blastdb.pl') + ' --decompress %s --passive --timeout 10' % self.db

            try:
                p = subprocess.run(arg, shell=True, stdout=subprocess.PIPE, cwd=self.dbpath)
                self.log.debug(p.stdout.decode('utf-8').strip())
            except subprocess.CalledProcessError as e:
                self.log.exception('BLAST+ db update failed, returned code %d' % e.returncode)

        # check database
        self.log.debug('checking %sBLAST db %s in %s'
                       % ('user-supplied ' if not self.update else '', self.db, self.dbpath))
        arg = '%s -db %s -info' % (path.join(self.blast_dir, 'blastdbcmd'),
                                   path.join(self.dbpath, self.db))
        try:
            res = subprocess.check_output(arg, shell=True).decode('utf-8')
            if 'error' in res:
                raise ValueError
        except (subprocess.CalledProcessError, ValueError) as e:
            self.log.error('Error in BLAST+ db')
            os._exit(1)  # this is no ordinary exit; it kills zombies, too!
        self.log.info('BLAST+ db was ok')

        # run BLAST+
        self.log.debug('BLASTing locally ...')
        arg = '%s -db %s -query %s -num_threads 3 -max_target_seqs 5 -outfmt 5 -out %s' \
              % (path.join(self.blast_dir, 'blastn'),
                 path.join(self.dbpath, self.db),
                 self.FASTA, self.XML)

        start = time()
        self.log.debug(arg)
        try:
            subprocess.run(arg, shell=True, check=True)
            self.log.info('finished BLAST+ in %.2f sec' % (time() - start))
        except subprocess.CalledProcessError as e:
            self.log.exception('BLAST failed, returned %d\n%s'
                               % (e.returncode, e.output.decode('utf-8') if e.output is not None else ''))
            return None

        # read in metadata tsv
        df = pandas.read_csv(self.TSV, sep='\t', dtype={'id': str})
        df.set_index('id', inplace=True)
        df['BLAST_species'] = ''
        df['pid'] = nan

        # parse .xml and write best hits info to metadata
        missing_seqs = list()
        for entry in SearchIO.parse(self.XML, 'blast-xml'):
            try:
                res = self._parse(entry)
                df.at[entry.id, 'BLAST_species'] = res[0]
                df.at[entry.id, 'pid'] = res[1]
            except ValueError:
                # no hits found -> save seq_id for online NCBI BLAST
                missing_seqs.append(entry.id)

        # write to TSV
        df.to_csv(self.TSV, sep='\t', na_rep='', header=True, index=True)
        self.log.debug('wrote updated metadata to %s' % self.TSV)

        return missing_seqs

    def _blast_remote(self, missing_seqs):
        """
        Write seqs missing from local BLAST+ db to .FASTA and search in NCBI nucleotides.
        Updates non-seq data .csv again.
        :return:
        """
        self.log.warning('online BLAST: %d seq%s missing from %s'
                         % (len(missing_seqs), 's' if len(missing_seqs) > 1 else '', self.db))

        # run remote NCBI BLAST via BLAST+
        self.log.debug('BLASTing online ...')
        arg = '%s -db %s -query %s -remote -max_target_seqs 5 -outfmt 5 -out %s' \
              % (path.join(self.blast_dir, 'blastn'), self.remote_db, self.missing_fasta, self.www_XML)

        start = time()
        self.log.debug(arg)
        try:
            subprocess.run(arg, shell=True, check=True, timeout=self.timeout)
            self.log.info('finished remote BLAST in %.2f sec' % (time() - start))
            self._parse_remote_result()
        except subprocess.CalledProcessError as e:
            self.log.exception('remote BLAST failed, returned %d\n%s'
                               % (e.returncode, e.output.decode('utf-8') if e.output is not None else ''))
        except subprocess.TimeoutExpired as e:
            self.log.error('remote BLAST timed out')

    def _parse(self, xml_entry):
        """Parses genus and species from a 'Hit' in the XML entry. Raises a ValueError if no hits."""
        try:
            _def = xml_entry.hits[0].description.split(',')[0].strip().split(' ')
            genus, species = _def[0], _def[1]

            if species == 'cf.':
                species += _def[2]
                self.log.warning('found a cf. but that\'s fine')

            # parse % identity from high-scoring pair
            _hsp = xml_entry.hsps[0]
            pid = _hsp.ident_num * 100 / (_hsp.aln_span + _hsp.gap_num)

            return ' '.join([genus, species]), pid
        except IndexError:
            raise ValueError

    def _parse_remote_result(self):
        """
        Parses an XML from a previous online BLAST run
        :return:
        """
        # read in metadata tsv
        df = pandas.read_csv(self.TSV, sep='\t', dtype={'id': str})
        df.set_index('id', inplace=True)

        # parse .xml and write best hits info to metadata
        with open(self.bad_seqs, 'a') as fh:
            for entry in SearchIO.parse(self.www_XML, 'blast-xml'):
                try:
                    res = self._parse(entry)
                    df.at[entry.id, 'BLAST_species'] = res[0]
                    df.at[entry.id, 'pid'] = res[1]
                except ValueError:
                    fh.write('%s\t%s\t%s\t%s\tno hit in NCBI BLAST nucleotide database\n'
                             % (df.at[entry.id, 'file'], entry.id,
                                df.at[entry.id, 'box'], self.gene))
                    self.log.error('%s no hit in NCBI nucleotide db' % entry.id)

        # write to TSV
        df.to_csv(self.TSV, sep='\t', na_rep='', header=True, index=True)
        self.log.debug('wrote newly updated metadata to %s' % self.TSV)
