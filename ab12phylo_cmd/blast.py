# 2021 Leo Kaindl

"""
Runs a BLAST+ search on a local database and a remote NCBI BLAST search in the nucleotides
database for seqs missing from it.
Result :code:`XML` files are parsed, and dictionary of non-sequence data is updated and
written to file. Processes are run in a separate thread due to the Python GIL.
"""

import logging
import multiprocessing
import os
import shutil
import sys
import urllib.error
from math import ceil
from os import path
from pathlib import Path
from subprocess import run, check_output, PIPE, CalledProcessError
from time import time, sleep

import pandas as pd
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from numpy import nan

from ab12phylo_cmd import cli, phylo
from ab12phylo_cmd.filter import chmod_x


def _add_xml(*args):
    """Entry point for adding BLAST results"""
    logger = logging.getLogger(__name__)
    # if no argument -> default to cwd
    args = sys.argv[1:]
    if len(args) == 0:
        args = [os.getcwd()]
    args = ['--add_xml'] + args
    # parse to get defaults
    namespace = cli.parser(args, visualize=True).args

    # find and read log
    try:
        log = open(namespace.log[:-5] + '1.log', 'r').read()
    except FileNotFoundError:
        try:
            log = open(namespace.log[:-7] + '.log', 'r').read()
        except FileNotFoundError:
            logger.error('no log file found')
            print('no log file found', file=sys.stderr)
            exit(1)

    # extract genes
    start = log.find('--GENES--') + 10
    end = log[start:start + 800].find('\n')
    namespace.genes = log[start:start + end].split('::')
    logger.debug('read genes from log: %s' % ':'.join(namespace.genes))
    print('read genes from log: %s' % ':'.join(namespace.genes))

    blaster = blast_build(namespace, None)
    print('start reading ...')
    blaster.start()
    # reading an XML requires no waiting
    blaster.join()
    # visualize MSA for inspection
    print('write HTML of MSA ...')
    phylo.mview_msa(namespace)
    print('BYE!')


class blast_build(multiprocessing.Process):
    """
    -none: do nothing
    -xml: read in result files and nothing else
    -local: do not run remote BLAST
    -remote: do not run local BLAST
    If a local search is run:
    Looks for pre-installed BLAST+ 2.9.0 or newer, updates or downloads a database if not
    specified via path, else checks integrity. Runs online Bio.BLAST.NCBIWWW.qblast unless
    prohibited, otherwise checks for earlier result.
    Writes information from best hit to metadata.tsv.
    """

    def __init__(self, _args, reader):
        multiprocessing.Process.__init__(self)
        self._stop_event = multiprocessing.Event()
        self.log = logging.getLogger(__name__)

        # modes / variables
        self.no_BLAST = _args.no_BLAST
        if self.no_BLAST is True:
            return
        self.no_remote = _args.no_remote
        self.no_local = _args.no_local
        self.gene = _args.genes[0]
        self.db = _args.db
        self.remote_db = _args.remote_db
        self.timeout = _args.timeout
        self.gui = 'gui' in _args
        self.df = pd.DataFrame()
        self.update = False

        # Files / Paths
        self.FASTA = path.join(_args.dir, self.gene, self.gene + '.fasta')
        self.XML = _args.xml
        self.www_XML = _args.www_xml
        self.read_XML = _args.BLAST_xml
        self.TSV = _args.tsv
        self.bad_seqs = _args.bad_seqs
        self.missing_fasta = _args.missing_fasta
        self.dbpath = _args.dbpath
        self.blast_dir = None
        os.makedirs(path.join(_args.dir, 'BLAST'), exist_ok=True)
        self.cfg = _args.cfg

        # data in memory
        try:
            self.seqdata = reader.seqdata
            self.metadata = reader.metadata
        except AttributeError as ex:
            self.log.info('passed empty reader')

    def run(self):
        if self.no_BLAST is True:
            self.log.info('--SKIPPING BLAST--')  # leave it alone!
            return
        while not self._stop_event.wait(.1):
            if self.read_XML:
                self.log.info('reading BLAST xml')
                try:
                    not_found = self._parse_remote_result(self.read_XML)
                    self._write_missing(not_found)
                except Exception as ex:
                    self.log.exception(ex)
                return

            if not self.no_local:
                # run local BLAST
                not_found = self._blast_local()
                # if local BLAST failed, continue without annotation
                if not_found is None:
                    return
            else:
                self.log.debug('skip local BLAST')
                # get all IDs for a remote BLAST
                not_found = list(self.seqdata[self.gene].keys())

            self._write_missing(not_found)

            if self.no_remote is True:
                self.log.info('skip online BLAST for missing seqs')
                # try reading old result
                if os.path.isfile(self.www_XML):
                    try:
                        self.log.warning('reading XML result of previous online BLAST')
                        self._parse_remote_result(self.www_XML)
                    except Exception as ex:
                        self.log.error(ex)
                return

            # online BLAST for seqs missing from local db
            self._blast_remote(not_found)
            return

    def stop(self):
        self.terminate()
        self._stop_event.set()
        self.log.warning('killed BLAST')

    def _blast_local(self):
        """
        Optionally check BLAST+ installation, update a local BLAST+ database, check integrity.
        Run a BLAST+ search and finally write updated non-seq data to a .csv file.
        :return: missing sequences as list of sample IDs
        """
        while not self._stop_event.wait(.1):
            # check BLAST+
            try:
                try:
                    binary = self.cfg['blastn']
                except KeyError:
                    binary = shutil.which('blastn')
                if binary is None:
                    raise ValueError('BLAST+ not installed (not on the $PATH)')
                chmod_x(binary)
                output = check_output(f'{binary} -version', shell=True).decode('utf-8')
                version = [int(i) for i in output.split(',')[0].split(' ')[-1].split('.')]
                if version[0] == 2 and version[1] >= 9 or version[0] > 2:
                    self.blast_dir = path.dirname(binary)
                    blast_dir = Path(self.blast_dir)
                else:
                    raise ValueError(f'BLAST+ in {self.blast_dir} is outdated')
            except (CalledProcessError, ValueError) as e:
                self.log.error(e)
                if self.gui:
                    raise e
                else:
                    exit(e)

            # if BLAST db path was provided, use and don't update
            if self.dbpath is not None:
                self.update = False
            else:
                # save BLAST db in root directory
                self.dbpath = path.join(path.dirname(path.dirname(__file__)), 'blastdb', self.db)
                self.update = True
                os.makedirs(self.dbpath, exist_ok=True)

            # update database
            if self.update:
                self.log.debug('downloading or updating BLAST+ db ...')
                binary = blast_dir / 'update_blastdb.pl'
                chmod_x(binary)
                arg = f'{binary} --decompress {self.db} --passive --timeout {self.timeout}'
                if sys.platform in ['win32', 'cygwin']:
                    arg = f"{shutil.which('perl')} {arg}"
                try:
                    self.log.debug(arg)
                    if sys.platform in ['win32', 'cygwin']:
                        p = run(arg, shell=True, stdout=PIPE, cwd=self.dbpath, creationflags=0x8000000)
                    else:
                        p = run(arg, shell=True, stdout=PIPE, cwd=self.dbpath)
                    self.log.debug(p.stdout.decode('utf-8').strip())
                except CalledProcessError as e:
                    self.log.exception('BLAST+ db update failed, returned code %d' % e.returncode)

            # check database
            self.log.debug('checking %sBLAST db %s in %s'
                           % ('user-supplied ' if not self.update else '', self.db, self.dbpath))
            binary = blast_dir / 'blastdbcmd'
            chmod_x(binary)
            db = path.join(self.dbpath, self.db)
            arg = f'{binary} -db "{db}" -info'
            try:
                res = check_output(arg, shell=True).decode('utf-8')
                if 'error' in res:
                    raise ValueError
            except (CalledProcessError, ValueError) as e:
                self.log.error('Error in BLAST+ db')
                os._exit(1)  # this is no ordinary exit; it kills zombies, too!
            self.log.debug('BLAST+ db was ok')

            # run BLAST+
            self.log.debug('BLASTing locally ...')
            binary = blast_dir / 'blastn'
            chmod_x(binary)
            arg = f'{binary} -db "{db}" -query "{self.FASTA}" -num_threads 3 ' \
                  f'-max_target_seqs 10 -outfmt 5 -out "{self.XML}"'
            start = time()
            self.log.debug(arg)
            try:
                if sys.platform in ['win32', 'cygwin']:
                    run(arg, shell=True, check=True, creationflags=0x8000000)
                else:
                    run(arg, shell=True, check=True)
                self.log.info(f'finished BLAST+ in {(time() - start):.2f} sec')
            except CalledProcessError as e:
                self.log.exception('BLAST failed, returned %d\n%s'
                                   % (e.returncode, e.output.decode('utf-8') if e.output is not None else ''))
                return None

            # read in metadata tsv
            self.df = pd.read_csv(self.TSV, sep='\t', dtype={'id': str})
            self.df.set_index(['id', 'gene'], inplace=True)
            self.df['BLAST_species'] = ''
            self.df['pid'] = nan
            self.df['extra_species'] = ''

            # parse .xml and write best hits info to metadata
            missing_seqs = list()
            for entry in SearchIO.parse(self.XML, 'blast-xml'):
                if not self._parse(entry):
                    # no hits found -> save seq_id for online NCBI BLAST
                    missing_seqs.append(entry.id)

            # write to TSV
            self.df.to_csv(self.TSV, sep='\t', na_rep='', header=True, index=True)
            self.log.debug('wrote updated metadata to %s' % self.TSV)

            return missing_seqs  # regular

        return []  # on thread kill

    def _blast_remote(self, missing_seqs):
        """
        Write seqs missing from local BLAST+ db to .FASTA and search in NCBI nucleotides.
        Updates non-seq data .csv again.
        :return:
        """
        while not self._stop_event.wait(.1):
            self.log.warning('online BLAST: %d seq%s missing from %s'
                             % (len(missing_seqs), 's' if len(missing_seqs) > 1 else '', self.db))

            # run remote NCBI BLAST via Biopython
            runs = ceil(len(missing_seqs) / 10)
            self.log.debug('BLASTing online in %d runs.' % runs)
            # prep file paths
            self.www_XML = [self.www_XML[:-4] + '_%d.xml' % run for run in range(1, runs + 1)]

            sleep_time = 11 if runs > 1 else 0
            for run in range(runs):
                start = time()
                # create query search string
                records = '\n'.join([self.seqdata[self.gene][seq_id].format('fasta')
                                     for seq_id in missing_seqs[run * 10: run * 10 + 10]])
                try:
                    # raise urllib.error.URLError('testing')
                    handle = NCBIWWW.qblast(program='blastn', database=self.remote_db,
                                            sequence=records, format_type='XML', hitlist_size=10)
                    # save as .xml
                    with open(self.www_XML[run], 'w') as fh:
                        fh.write(handle.read())
                    self.log.info('remote BLAST %d:%d complete after %.2f sec' % (run + 1, runs, time() - start))
                    self._parse_remote_result([self.www_XML[run]])
                    sleep(sleep_time)
                except urllib.error.URLError as offline:
                    self.log.warning('online BLAST aborted, no internet connection')
                    return
                except Exception as e:
                    self.log.exception(e)
                    sleep(sleep_time)

            self.log.info('finished remote BLAST')
            return

    def _parse(self, xml_entry):
        """
        Parse genus, species, pid, share of hits and other species from XML entry,
        assuming the calculated pid will only fall.
        :param xml_entry:
        :return: True if successful, False if no hits or Error.
        """
        if not xml_entry:
            return False
        try:
            # get row for selected gene
            row = self.df.loc[xml_entry.id, self.gene]
        except KeyError:
            try:
                rows = self.df.loc[xml_entry.id]
                self.gene = [rows.index[0], self.gene]
                row = rows.iloc[0]
            except KeyError:
                return False
        try:
            sp_counts, sp_pids = dict(), dict()
            last_pid = 0

            for hit in xml_entry:
                _def = hit.description.split(',')[0].strip().split(' ')
                genus, species = _def[0], _def[1]

                if species == 'cf.':
                    species += _def[2]
                    self.log.warning('found a cf. but that\'s fine')
                species = f'{genus} {species}'

                # parse % identity from high-scoring pair
                _hsp = hit.hsps[0]
                pid = _hsp.ident_num * 100 / (_hsp.aln_span + _hsp.gap_num)

                sp_pids[species] = pid
                sp_counts[species] = sp_counts.get(species, 0) + 1

                if pid != last_pid and len(set(sp_pids.values())) >= 2:
                    break  # as soon as there are two different pids, there must be two different species
                last_pid = pid

            finders = list()
            if len(sp_pids) == 2:
                # easiest case: found two species with two different pids
                four = [f for fs in [[k, v] for k, v in sp_pids.items()] for f in fs]
                if four[3] > four[1]:
                    finders.append([four.pop(3), four.pop(2)])
                while four:
                    finders.append([four.pop(1), four.pop(0)])
            else:
                # iterate over observed species along decreasing pid
                for best_pid in sorted(set(sp_pids.values()), reverse=True):
                    sp_best_pid = [sp for sp, pid in sp_pids.items() if pid == best_pid]
                    if len(sp_best_pid) == 1:
                        # best pid occurs exactly once
                        finders.append([best_pid, sp_best_pid[0]])
                    else:
                        # several species have best pid
                        sp_counts = {sp: count for sp, count in sp_counts.items() if sp in sp_best_pid}
                        most_freq = {sp for sp, c in sp_counts.items() if c == max(sp_counts.values())}
                        finders.append([best_pid, ' / '.join(most_freq)])
            if not finders:
                return False

            # write to dataframe
            for i, j in enumerate(['pid', 'BLAST_species']):
                row.at[j] = finders[0][i]
            if len(finders) > 1:
                row.at['extra_species'] = '%.2f: %s' % (finders[1][0], finders[1][1])
            if type(self.gene) == list:
                gene, self.gene = self.gene
                self.df.loc[xml_entry.id, gene] = row
            else:
                self.df.loc[xml_entry.id, self.gene] = row

            return True
        except Exception as ex:
            self.log.exception(ex)
            return False

    def _parse_remote_result(self, in_xmls):
        """
        Parses a list of XML files from previous online BLAST run
        :return:
        """
        # file check
        xmls = list()
        for xml in in_xmls:
            if path.isfile(xml):
                xmls.append(xml)
            else:
                self.log.warning('XML path invalid: %s' % xml)
        if not xmls:
            self.log.warning('no valid paths.')
            return
        # read in metadata tsv
        self.df = pd.read_csv(self.TSV, sep='\t', dtype={'id': str, 'BLAST_species': str})
        self.df.set_index(['id', 'gene'], inplace=True)
        if 'BLAST_species' not in self.df.columns:
            self.df['BLAST_species'] = ''
            self.df['pid'] = nan
            self.df['extra_species'] = ''

        errors = False
        # parse .xml and write best hits info to metadata
        with open(self.bad_seqs, 'a') as fh:
            for xml in xmls:
                for entry in SearchIO.parse(xml, 'blast-xml'):
                    if not self._parse(entry):
                        try:
                            if entry.id in self.df.index:
                                file, box = self.df.loc[entry.id, self.gene][['file', 'box']]
                                fh.write(f'{file}\t{entry.id}\t{box}\t{self.gene}\t'
                                         f'no hit in NCBI {self.remote_db} database\n')
                                self.log.error(f'{entry.id} no hit in NCBI {self.remote_db} db')
                            else:
                                self.log.error(f'no entry {entry.id} in metadata TSV')
                        except KeyError as ke:
                            errors = True
                            self.log.error(ke)
        if errors:
            self.log.warning('printing DataFrame for check')
            self.log.warning(self.df)

        # write to TSV
        self.df.to_csv(self.TSV, sep='\t', na_rep='', header=True, index=True)
        self.log.debug('wrote newly updated metadata to %s' % self.TSV)

        not_found = \
            set(self.df[self.df.BLAST_species == ''].index.get_level_values(0)) - \
            set(self.df[self.df.BLAST_species != ''].index.get_level_values(0))
        return not_found

    def _write_missing(self, not_found):
        """Save missing sequences to file"""
        if not not_found:
            self.log.info('all sequences found')
            return
        with open(self.missing_fasta, 'w') as fh:
            for seq_id in not_found:
                for gene in [self.gene] + list(self.seqdata.keys() - {self.gene}):
                    r = self.seqdata[gene].get(seq_id, False)
                    if r:
                        fh.write(r.format('fasta') + '\n')
                        break
