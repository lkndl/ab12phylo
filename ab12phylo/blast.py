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
import urllib.error
from os import path
from time import time

import pandas
from Bio import SearchIO
from Bio.Blast import NCBIWWW


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
        elif self.no_remote is True:
            self.log.warning('skip online BLAST for missing seqs')

            # try reading old result
            if os.path.isfile(self.www_XML):
                try:
                    self.log.warning('reading XML result of previous NCBI BLAST')
                    self._parse_remote_result()
                except Exception as ex:
                    self.log.error(ex)
            return

        self._blast_remote(not_found)
        self._parse_remote_result()

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

        # parse .xml and write best hits info to metadata
        missing_seqs = list()
        for entry in SearchIO.parse(self.XML, 'blast-xml'):
            try:
                self.metadata[self.gene][entry.id].update(self._parse(entry))
            except ValueError:
                # no hits found -> save seq_id for online NCBI BLAST
                missing_seqs.append(entry.id)

        # make a pandas DataFrame from metadata
        df = pandas.concat({gene: pandas.DataFrame.from_dict(records, orient='index')
                            for gene, records in self.metadata.items()})
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

        # create query seq string
        records = '\n'.join([self.seqdata[self.gene][seq_id].format('fasta') for seq_id in missing_seqs])

        # save to file
        with open(self.missing_fasta, 'w') as fh:
            fh.write(records)
            self.log.debug('sequences not found in %s saved in %s' % (self.db, self.missing_fasta))

        # send query
        self.log.info('BLASTing online ...')
        start = time()
        self.log.info('sequences not found in %s: %s' % (self.db, ' '.join(missing_seqs)))
        try:
            handle = NCBIWWW.qblast(program='blastn', database='nt',
                                    sequence=records, format_type='XML', hitlist_size=10)
            # save as .xml
            self.log.info('finished NCBI BLAST in %.2f sec' % (time() - start))
            with open(self.www_XML, 'w') as xml:
                xml.write(handle.read())
        except urllib.error.URLError as offline:
            self.log.warning('online BLAST aborted, no internet connection')
            return
        except Exception as e:
            self.log.exception(e)
            return

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

            return {'BLAST_species': ' '.join([genus, species]), 'pid': pid}
        except IndexError:
            raise ValueError

    def _parse_remote_result(self):
        """
        Parses an XML from a previous online BLAST run
        :return:
        """
        # parse .xml and write best hits info to metadata
        with open(self.bad_seqs, 'a') as fh:
            for entry in SearchIO.parse(self.www_XML, 'blast-xml'):
                try:
                    self.metadata[self.gene][entry.id].update(self._parse(entry))
                except ValueError:
                    fh.write('%s\t%s\t%s\t%s\tno hit in NCBI BLAST nucleotide database\n'
                             % (self.metadata[self.gene][entry.id]['file'], entry.id,
                                self.metadata[self.gene][entry.id]['box'], self.gene))
                    self.log.error('%s no hit in NCBI nucleotide db' % entry.id)

        # make a pandas DataFrame from metadata
        df = pandas.concat({gene: pandas.DataFrame.from_dict(records, orient='index')
                            for gene, records in self.metadata.items()})
        # write to TSV
        df.to_csv(self.TSV, sep='\t', na_rep='', header=True, index=True)
        self.log.debug('wrote newly updated metadata to %s' % self.TSV)
