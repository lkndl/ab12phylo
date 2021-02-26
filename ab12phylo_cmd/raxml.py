# 2021 Leo Kaindl

"""
Builds phylogenetic trees using raxml-ng, parallelized in threads.
"""

import logging
import os
import random
import re
import shutil
import string
import subprocess
import threading
from os import path

import pandas


class raxml_build:
    """
    Builds phylogenetic trees using raxml-ng.
    https://doi.org/10.1093/bioinformatics/btz305
    """

    def __init__(self, args):
        self.log = logging.getLogger(__name__)
        self.args = args
        self.runs = self.cpus = self._prefixes = self._best_tree = None
        self.msa = self.args.msa

        # seed RNG
        random.seed(args.seed)

        # make a directory for raxml
        self._dir = path.join(args.dir, 'raxml')
        os.makedirs(self._dir, exist_ok=True)

        # look for raxml-ng binary
        self._binary = shutil.which('raxml-ng')
        local = True
        if self._binary is None:
            local = False
            self._binary = path.join(path.abspath(path.dirname(__file__)),
                                     'tools', 'raxml-ng-static', 'raxml-ng').replace(' ', '\ ')
        self.log.debug('use %s raxml-ng' % ('pre-installed' if local else 'out-of-the-box'))

    def run(self):
        """
        Calls MSA check. Then creates ML search threads and keeps track of used prefixes.
        Finds best tree topology from logs, then runs Bootstrapping threads.
        Computes support values for best topology from collated bootstrap trees.
        Finally copies supported tree to parent directory.
        :return:
        """
        # check+parse msa, determine number of parallel raxml-ng instances and cpus for each
        self.runs, self.cpus, self.msa = self._check_msa

        # keep track of finished prefixes:
        self._prefixes = []

        # keep track of threads:
        active_threads = []

        self.log.debug('starting ML searches')
        # calculate ML tree searches share
        ml_searches = [mltype // self.runs + 1 for mltype in self.args.start_trees]

        # create+start ML search threads
        for j in range(self.runs):
            prefix = 'ml_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=4))
            self._prefixes.append(prefix)
            seed = random.randint(0, max(1000, self.args.bootstrap))
            thread = raxml_thread((prefix, self.log, self._binary, self.msa, ml_searches[0],
                                   ml_searches[1], self.args.evomodel, self.args.min_dist, seed,
                                   path.join(self._dir, prefix), self.cpus), mode='infer_topology')
            active_threads.append(thread)
            thread.start()
        self.log.info('active Threads: %s' % ' '.join(self._prefixes))

        # finish threads
        for thread in active_threads:
            thread.join()
            self.log.info('finished %s' % thread.prefix)
        self.log.debug('all ML searches done')

        # find best topology: parse logs for best LogLikelihood
        best_prefix, best_score = '', float('-inf')
        for prefix in self._prefixes:
            with open(path.join(self._dir, prefix + '.raxml.log'), 'r') as log:
                text = log.read()
                pos = text.find('Final LogLikelihood:')
                score = float(text[pos:pos + 40].split('\n')[0].split(' ')[-1])
                if score > best_score:
                    best_prefix, best_score = prefix, score
        self.log.info('best topology %s %.2f' % (best_prefix, best_score))

        # set to file path
        self._best_tree = path.join(self._dir, best_prefix + '.raxml.bestTree')

        self.log.debug('starting Bootstrapping')
        self._prefixes = []
        active_threads = []
        # calculate bootstraps share
        trees = self.args.bootstrap // self.runs + 1

        # create+start Bootstrapping threads
        for j in range(self.runs):
            prefix = 'bs_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=4))
            self._prefixes.append(prefix)
            seed = random.randint(0, max(1000, self.args.bootstrap))
            thread = raxml_thread((prefix, self.log, self._binary, self.msa, self._best_tree,
                                   trees, self.args.evomodel, self.args.min_dist, seed,
                                   path.join(self._dir, prefix), self.cpus), mode='bootstrap')
            active_threads.append(thread)
            thread.start()
        self.log.info('active BS Threads: %s' % ' '.join(self._prefixes))

        # finish threads
        for thread in active_threads:
            thread.join()
            self.log.info('finished %s' % thread.prefix)

        self.log.debug('computing supports')
        # concat all bootstrap trees
        all_bs_trees = path.join(self._dir, 'all_bs_trees.nwk')
        with open(all_bs_trees, 'w') as fh:
            for prefix in self._prefixes:
                try:
                    with open(path.join(self._dir, prefix + '.raxml.bootstraps'), 'r') as bootstraps:
                        fh.write(bootstraps.read() + '\n')
                except FileNotFoundError:
                    self.log.warning('results from %s missing' % prefix)
        self.log.debug('all bootstrap trees in %s' % all_bs_trees)

        run_cmd = '%s --support --tree "%s" --bs-trees "%s" --prefix "%s" ' \
                  '--threads %d --bs-metric fbp,tbe --redo' \
                  % (self._binary, self._best_tree, all_bs_trees,
                     path.join(self._dir, '_sup'), self.cpus)

        _run_sp(run_cmd)

        # copy best trees with the right support values to parent directory
        shutil.move(path.join(self._dir, '_sup.raxml.supportFBP'), self.args.final_tree + '_FBP.nwk')
        shutil.move(path.join(self._dir, '_sup.raxml.supportTBE'), self.args.final_tree + '_TBE.nwk')
        self.log.info('final trees %s_[FBP,TBE].nwk' % self.args.final_tree)

    @property
    def _check_msa(self):
        """
        Checks if the MSA can be understood by raxml-ng,
        parses it to binary and estimates effective number of threads.
        If replacing is off, the original MSA will be used.
        :return: number of parallel instances possible on this machine,
        number of CPUs recommended by raxml-ng for each one, and new path to MSA.
        """
        # prep paths for checking, reducing and parsing MSA
        paths = [path.join(self._dir, prefix) for prefix in ['_chk', '_red']]

        # check msa
        arg = '%s --msa "%s" --check --model %s --prefix "%s"' \
              % (self._binary, self.msa, self.args.evomodel, paths[0])

        res = _run_sp(arg)

        keepers = set()
        if 'replace' in self.args and self.args.replace:
            # replace mode: replace sets of identical sequences by one representative. off by default
            # parse MSA check result for success and identical sequences
            pattern = re.compile(r'quences (.+?) and (.+?) are exactly identi')
            goners = 0

            gene = self.args.genes[0]
            df = pandas.read_csv(self.args.tsv, sep='\t', dtype={'id': str})
            df.set_index('id', inplace=True)
            df['replaces'] = ''
            df['replaced_by'] = ''

            # make temporary sub_df
            df1 = df[df.gene == gene]

            for keeper, goner in re.findall(pattern, res):
                goners += 1
                keepers.add(keeper)
                if df1.at[keeper, 'replaces'] == '':
                    df1.at[keeper, 'replaces'] = goner
                else:
                    df1.at[keeper, 'replaces'] = df1.at[keeper, 'replaces'] + ', ' + goner
                df1.at[goner, 'replaced_by'] = keeper

            # re-merge data frames
            df[df.gene == gene] = df1
            # swap columns
            cols = df.columns.to_list()

            # write to TSV
            df.to_csv(self.args.tsv, sep='\t', na_rep='', header=True, index=True)
            if goners > 0:
                self.log.warning('Identical sequences were found, %d replacements occurred' % goners)

        # look if MSA check was ok
        if 'successfully' not in res:
            self.log.error('error in MSA check, see %s' % paths[0] + '.raxml.log')
            exit(1)
        self.log.info('MSA check ok')

        # if any duplicates were found, a reduced alignment was written. ignore if not in replace mode
        msa = paths[0] + '.raxml.reduced.phy' if len(keepers) > 0 else self.msa

        # now parse, and use reduced alignment if duplicates found. resulting .rba will be smaller
        arg = '%s --msa "%s" --parse --model %s --prefix "%s"' \
              % (self._binary, msa, self.args.evomodel, paths[1])

        res = _run_sp(arg)

        # find number of threads:
        pos = res.find('* Recommended number of threads')
        cpus = int(res[pos:pos + 60].split('\n')[0].split(' ')[-1])

        # how many parallel runs on how many cores?
        runs = max(1, os.cpu_count() // cpus)
        # potentially limit threads
        if 'max_threads' in self.args:
            runs = min(runs, self.args.max_threads)
        self.log.info('use %d cpu%s per raxml-ng. found %d physical cores. create %d instances'
                      % (cpus, '' if cpus == 1 else 's', os.cpu_count(), runs))

        if 'replace' in self.args and self.args.replace:
            return runs, cpus, path.join(paths[1] + '.raxml.rba')
        else:
            # if not in replace mode, just return the unparsed alignment
            return runs, cpus, self.msa


def _run_sp(cmd):
    """
    Runs a command in a subprocess. Exceptions will be caught and logged, but cause exit.

    :param cmd: command that will be run
    :return: None
    """
    log = logging.getLogger(__name__)
    log.debug(cmd)
    try:
        return subprocess.check_output(cmd, shell=True).decode('utf-8')
    except subprocess.CalledProcessError as e:
        if e.returncode == -2:
            log.exception('RAxML user interrupt')
        else:
            log.exception('RAxML process failed, returned %s\n' % str(e.returncode)
                          + e.output.decode('utf-8') if e.output is not None else '')
        os._exit(1)  # this is no ordinary exit; it kills zombies, too!


class raxml_thread(threading.Thread):
    """A thread that creates raxml-ng trees"""

    def __init__(self, args, mode):
        threading.Thread.__init__(self)
        self.prefix = args[0]
        self.log = args[1]
        if mode == 'infer_topology':
            self.run_cmd = '%s --msa "%s" --tree rand{%d},pars{%d} --model %s --blmin %s ' \
                           '--seed %d --prefix "%s" --threads %d --redo' % args[2:]
        elif mode == 'bootstrap':
            self.run_cmd = '%s --bootstrap --msa "%s" --tree "%s" --bs-trees autoMRE{%d} ' \
                           '--model "%s" --blmin %s  ' \
                           '--seed %d --prefix "%s" --threads %d --redo' % args[2:]
        else:
            self.log.error('wrong raxml_thread mode')
            os._exit(1)  # this is no ordinary exit; it kills zombies, too!

    def run(self):
        try:
            subprocess.check_output(self.run_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            self.log.error('RAxML thread %s failed' % self.prefix)
