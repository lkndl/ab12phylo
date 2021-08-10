# 2021 Leo Kaindl

"""
This module builds a multiple sequence alignment for each gene and concats
them into one :code:`FASTA` file as input for the :class:`raxml` module.
"""

import logging
import os
import random
import shutil
import subprocess
import sys
from os import path
from pathlib import Path
from time import time, sleep

from Bio import SeqIO

from ab12phylo_cmd.filter import chmod_x


class msa_build:
    """Builds a Multiple Sequence Alignment"""

    def __init__(self, args, seq_counts, no_run=False):
        self.log = logging.getLogger(__name__)
        self.dir = args.dir
        self.genes = args.genes
        self.algo = args.msa_algo
        self.email = args.user
        self.msa = args.msa
        self.sep = args.sep
        self.missing_samples = args.missing_samples
        self.tools_path = path.join(path.abspath(path.dirname(__file__)), 'tools')

        # look for pre-installed version of selected algorithm
        self.binary = shutil.which(args.msa_algo)

        if no_run:
            return

        # build MSAs
        for gene in self.genes:
            if self.binary is not None:
                self.build_local(gene)
            else:
                self.build_remote(gene)
                sleep(1)

            # trim MSAs using Gblocks
            self.trim_msa(gene, seq_counts[gene], args.gblocks)

        # concat MSAs
        self.concat_msa()

    def reset_paths(self, _dir, msa):
        self.dir = _dir
        self.msa = msa

    def build_local(self, gene, no_run=False, new_arg=False):
        """Build MSAs locally using a pre-installed binary"""
        log_file = path.join(self.dir, gene, self.algo + '.log')
        fasta = path.join(self.dir, gene, gene + '.fasta')
        raw_msa = path.join(self.dir, gene, gene + '_raw_msa.fasta')
        self.log.debug('preparing %s MSA run' % self.algo)

        if new_arg:
            arg = new_arg

        elif self.algo == 'mafft':
            arg = '%s --thread %d --auto  "%s" > "%s"' \
                  % (self.binary, os.cpu_count(), fasta, raw_msa)

        elif self.algo == 'clustalo':
            arg = '%s --in "%s" --out "%s" --outfmt fasta --threads %d --force --verbose --auto' \
                  % (self.binary, fasta, raw_msa, os.cpu_count())
            if sys.platform in ['win32', 'cygwin']:
                arg += ' & exit /b 0'

        elif self.algo == 'muscle':
            arg = '%s -in "%s" -out "%s"' \
                  % (self.binary, fasta, raw_msa)

        elif self.algo == 't_coffee':
            arg = '%s -in "%s" -out "%s" -output fasta_aln -type dna ' \
                  % (self.binary, fasta, raw_msa)
        else:
            if no_run:
                arg = ''  # when we're only interested in the object and on windows
            else:
                assert False, 'invalid MSA algorithm selected: %s' % self.algo

        if no_run:
            return arg
        self._run(arg, log_file, 'pre-installed %s' % self.algo, no_exit=True)

    def build_remote(self, gene, no_run=False, new_arg=False):
        """Builds an MSA online using an EBI API client"""
        log_file = path.join(self.dir, gene, self.algo + '.log')
        fasta = path.join(self.dir, gene, gene + '.fasta')
        raw_msa = path.join(self.dir, gene, gene + '_raw_msa.fasta')

        self.log.warning(f'{"preparing" if no_run else "running"} '
                         f'online {self.algo.upper()}')

        if new_arg:
            arg = new_arg
        else:
            if self.algo == 't_coffee':
                self.algo = 'tcoffee'

            # create base call
            arg = '"%s" "%s" --email %s --outfile "%s" --sequence "%s" ' \
                  % (sys.executable, path.join(self.tools_path, 'MSA_clients', self.algo + '.py'),
                     self.email, path.join(self.dir, gene, 'msa'), fasta)

            # adapt for specific algorithm
            if self.algo == 'mafft':
                arg += '--stype dna'
            elif self.algo == 'clustalo':
                arg += '--stype dna --outfmt fa'
            elif self.algo == 'muscle':
                arg += '--format fasta'
            elif self.algo == 'tcoffee':
                arg += '--stype dna --format fasta_aln'

        if no_run:
            return arg

        # build an MSA for each gene
        self._run(arg, log_file, 'online %s' % self.algo)
        try:
            shutil.move(path.join(self.dir, gene, 'msa.aln-fasta.fasta'), raw_msa)
        except FileNotFoundError:
            shutil.move(path.join(self.dir, gene, '%s_raw_msa.fasta' % gene), raw_msa)

    def trim_msa(self, gene, seq_count, gblocks_mode):
        """
        Trims an MSA using Gblocks, using a pre-installed or deployed version.

        :param gene: this helps find the right files to trim
        :param seq_count: for computing settings
        :param gblocks_mode: can be ['skip', 'relaxed', 'balanced', 'semi_strict', 'strict']
        :return:
        """
        log_file = path.join(self.dir, gene, 'gblocks.log')
        raw_msa = path.join(self.dir, gene, gene + '_raw_msa.fasta')

        if gblocks_mode == 'skip':
            shutil.copy(raw_msa, path.join(self.dir, gene, gene + '_msa.fasta'))
            self.log.info('skipped Gblocks trimming, only copied file')

        else:
            # look for local Gblocks
            binary = shutil.which('Gblocks')
            local = True
            if binary is None:
                # pick deployed Gblocks
                binary = path.join(self.tools_path, 'Gblocks_0.91b', 'Gblocks').replace(' ', '\ ')
                local = False
                if sys.platform in ['win32', 'cygwin']:
                    binary += '.exe'
            chmod_x(binary)

            # set Gblocks options
            b4 = 5
            if gblocks_mode == 'relaxed':
                # set the minimal permissible minimum number of identical
                # sequences per position to define a conserved position
                cons = seq_count // 2 + 1
                # minimal number for a conserved flanking position
                flank = cons
                # keep no, half or all gap positions
                gaps = ['n', 'h', 'a'][1]
                self.log.info('running relaxed Gblocks')

            elif gblocks_mode == 'balanced':
                cons = seq_count // 2 + 1
                flank = min(seq_count // 4 * 3 + 1, seq_count)
                gaps = ['n', 'h', 'a'][1]
                self.log.info('running balanced Gblocks')

            elif gblocks_mode == 'default':
                cons = seq_count // 2 + 1
                flank = min(int(seq_count * 0.85) + 1, seq_count)
                gaps = ['n', 'h', 'a'][0]
                b4 = 10
                self.log.info('running Gblocks at default settings')

            else:
                cons = int(seq_count * 0.9)
                flank = cons
                gaps = ['n', 'h', 'a'][0]
                self.log.info('running strict Gblocks')

            # create base call
            arg = '"%s" "%s" -t=d -b2=%d -b1=%d -b4=%d -b5=%s -e=".txt" -d=n -s=y -p=n' \
                  % (binary, Path(raw_msa).resolve(), flank, cons, b4, gaps)  # don't swap order!
            # force return code
            if sys.platform in ['win32', 'cygwin']:
                arg += ' & exit /b 0'
            else:
                arg += '; exit 0'
            # MARK the -t=d sets the mode to nucleotides ... adapt?
            self._run(arg, log_file, 'pre-installed Gblocks' if local else 'out-of-the-box Gblocks')
            shutil.move(raw_msa + '.txt', path.join(self.dir, gene, gene + '_msa.fasta'))

    def concat_msa(self, gui=False):
        """Reads all trimmed MSAs to memory, then iterates over samples, writes concatenated MSA"""
        self.log.debug('concatenating per-gene MSAs')
        # read in all MSAs using SeqIO
        all_records = {gene: {record.id: record.upper() for record in SeqIO.parse(
            path.join(self.dir, gene, gene + '_msa.fasta'), 'fasta')} for gene in self.genes}

        # get the length of the trimmed concat MSA
        msa_len = 0
        for gene in all_records.keys():
            msa_len += len(random.choice(list(all_records[gene].values())))

        missing_genes = {gene: list() for gene in self.genes}

        with open(self.msa, 'w') as msa:
            shared = 0
            # iterate over samples available for first gene
            _ids = gui if gui else set(all_records[self.genes[0]])
            for sample_id in _ids:
                # get SeqRecord for first gene
                record = all_records[self.genes[0]].pop(sample_id)

                skip = False
                # append other genes
                for gene in self.genes[1:]:
                    try:
                        record += self.sep  # to visually separate genes in the MSA
                        record += all_records[gene].pop(sample_id)
                    except KeyError:
                        missing_genes[gene].append(sample_id)
                        skip = True
                # write to file
                if not skip:
                    shared += 1
                    SeqIO.write(record, msa, 'fasta')
        if len(self.genes) > 1:
            if shared == 0:
                msg = 'No samples shared across all genes.'
                self.log.error(msg)
                raise ValueError(msg) if gui else exit(1)
            else:
                self.log.info('finished writing concat MSA with %d entries' % shared)
        else:
            if msa_len == 0:
                msg = 'No conserved sites found. Please try a more relaxed trimming mode!'
                self.log.error(msg)
                if not gui:
                    raise ValueError(msg) if gui else exit(1)
            self.log.info('copied MSA to result root')
        self.log.info('MSA shape: %dx%d' % (msa_len, shared))

        if gui:
            return msa_len, shared

        # any remaining samples were missing from first gene
        [missing_genes[self.genes[0]].extend(all_records[gene].keys()) for gene in self.genes[1:]]

        # write info about missing samples to file and log
        missing = {gene: ', '.join(set(samples)) for gene, samples
                   in missing_genes.items() if set(samples)}
        if missing:
            with open(self.missing_samples, 'w') as fh:
                fh.write('gene\tmissing samples\n')
                for k, v in missing.items():
                    fh.write('%s\t%s\n' % (k, v))
                    self.log.info('samples missing from %s: %s' % (k, v))

    def _run(self, arg, log_file, info, no_exit=False):
        """Runs a command and writes output to log-file"""
        self.log.debug(arg)
        start = time()
        with open(log_file, 'w') as log_handle:
            try:
                subprocess.run(arg, shell=True, check=True, stdout=log_handle, stderr=log_handle)
            except (OSError, subprocess.CalledProcessError) as e:
                if not no_exit:
                    self.log.exception('returned: ' + str(e.returncode)
                                       + e.output.decode('utf-8') if e.output is not None else '')
                    sys.exit(1)
                else:
                    raise e
        self.log.debug('%.2f seconds for %s' % (time() - start, info))
