# 2021 Leo Kaindl

"""
Builds phylogenetic trees using raxml-ng or iqtree2, parallelized in threads.
"""

import logging
import mmap
import os
import random
import re
import shutil
import string
import subprocess
import sys
import threading
from pathlib import Path

import pandas

from ab12phylo_cmd.filter import chmod_x


def ml_build(args):
    """This is just an intermediary to keep the code in main.py a tiny bit more concise."""
    if args.ml_tool == 'iqtree2' or sys.platform == 'win32':
        return iqtree_build(args)
    else:
        return raxml_build(args)


class iqtree_build:
    """
    Builds phylogenetic tree using iqtree-2.
    https://doi.org/10.1093/molbev/msu300
    """

    def __init__(self, args):
        self.log = logging.getLogger(__name__)
        self.args = args

        if self.args.ultrafast and self.args.bootstrap < 1000:
            self.log.warning(f'The selected number of iterations ({self.args.bootstrap}) '
                             f'is too low for ultrafast bootstrapping. '
                             f'Will run 1000 iterations.')
            self.args.bootstrap = 1000

        # make a directory for iqtree2
        self._dir = Path(args.dir) / 'iqtree'
        Path.mkdir(self._dir, exist_ok=True)

        # look for iqtree binary
        try:
            self._binary = self.args.cfg.get('iqtree2', shutil.which('iqtree2'))
        except KeyError:
            self._binary = shutil.which('iqtree2')
        if self._binary is None:
            self.log.error('IQ-Tree 2 not installed')
            os._exit(1)  # this is no ordinary exit; it kills zombies, too!
        # Make sure it's executable
        chmod_x(self._binary)

    def run(self):

        calls = list()
        if self.args.findmodel:
            calls.append('# find best model\n'
                         '"{iqtree2}" -s "{msa}" -m MF -mtree -mset raxml '
                         '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                         '--prefix "{mf_prefix}" ')

        calls.append('# infer ML tree\n'
                     '"{iqtree2}" -s "{msa}" -m "{evo_model}" '
                     '-ninit {start_trees} -ntop {start_trees} '
                     '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                     '--prefix "{ml_prefix}" ')

        if not self.args.ultrafast:
            calls.append('# non-parametric bootstrapping\n'
                         '"{iqtree2}" -s "{msa}" -m "{evo_model}" '
                         '-te "{ml_prefix}.treefile" -b {bootstraps} '
                         '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                         '--prefix "{bs_prefix}" ')
        else:
            calls.append('# ultrafast bootstrapping\n'
                         '"{iqtree2}" -s "{msa}" -m "{evo_model}" '
                         '-t "{ml_prefix}.treefile" -B {bootstraps} -wbtl '
                         '--seed {seed} -nt AUTO -ntmax {cpus} -redo '
                         '--prefix "{uf_prefix}" ')

        calls.append('# calc. TBE branch support\n'
                     '"{iqtree2}" -sup "{ml_prefix}.treefile" '
                     '-t "{boot_trees}" --tbe '
                     '-seed {seed} -nt AUTO -ntmax {cpus} -redo '
                     '--prefix "{sp_prefix}" ')

        # potentially limit threads
        cpus = os.cpu_count()
        if 'max_threads' in self.args:
            cpus = min(cpus, self.args.max_threads)

        # prep formatting args
        fargs = {'iqtree2': self._binary,
                 'msa': self.args.msa,
                 'seed': self.args.seed,
                 'evo_model': self.args.evomodel,
                 'cpus': cpus,
                 'mf_prefix': self._dir / 'mf',
                 'ml_prefix': self._dir / 'ml',
                 'bs_prefix': self._dir / 'bs',
                 'uf_prefix': self._dir / 'uf',
                 'sp_prefix': self._dir / 'sp',
                 'bootstraps': self.args.bootstrap,
                 'boot_trees': self._dir / ('uf.ufboot' if self.args.ultrafast else 'bs.boottrees'),
                 'start_trees': sum(self.args.start_trees)}

        for i, (desc, call) in enumerate([c.split('\n') for c in calls]):
            self.log.info(desc[2:])
            _run_sp(call.format(**fargs))

            if self.args.findmodel and i == 0:
                with open(f'{self._dir / "mf"}.iqtree') as fh:
                    s = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
                    pos = s.find(b'Best-fit model')
                    if pos != -1:
                        model = s[pos + 10:pos + 190].decode().split('\n')[0].split(':')[-1].strip()
                        msg = f'Extracted model: {model}'
                        self.log.info(msg)
                        fargs['evo_model'] = model
                    else:
                        self.log.error(f'Could not extract model from\n{self._dir / "mf"}.iqtree')
                        os._exit(1)  # this is no ordinary exit; it kills zombies, too!

        # copy trees with the right support values to parent directory
        shutil.copy(f'{self._dir / ("uf" if self.args.ultrafast else "bs")}.contree',
                    self.args.final_tree + '_FBP.nwk')
        shutil.copy(self._dir / 'sp.suptree', self.args.final_tree + '_TBE.nwk')
        self.log.info('final trees %s_[FBP,TBE].nwk' % self.args.final_tree)


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
        self._dir = Path(args.dir) / 'raxml'
        Path.mkdir(self._dir, exist_ok=True)

        # look for raxml-ng binary
        try:
            self._binary = self.args.cfg.get('raxml-ng', shutil.which('raxml-ng'))
        except KeyError:
            self._binary = shutil.which('raxml-ng')
        if self._binary is None:
            self.log.error('RAxML-NG not installed')
            os._exit(1)  # this is no ordinary exit; it kills zombies, too!
        # Make sure it's executable
        chmod_x(self._binary)

        if self.args.ultrafast or self.args.findmodel:
            self.log.warning('IQ-Tree parameters were passed for a RAxML-NG run, will be ignored.')

    def run(self):
        """
        Calls MSA check. Then creates ML search threads and keeps track of used prefixes.
        Finds best tree topology from logs, then runs Bootstrapping threads.
        Computes support values for best topology from collated bootstrap trees.
        Finally copies supported tree to parent directory.
        :return:
        """
        # check+parse msa, determine number of parallel raxml-ng instances and cpus for each
        self.runs, self.cpus, self.msa = self._check_msa()

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
                                   self._dir / prefix, self.cpus), mode='infer_topology')
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
            with open(f'{self._dir / prefix}.raxml.log', 'r') as log:
                text = log.read()
                pos = text.find('Final LogLikelihood:')
                score = float(text[pos:pos + 40].split('\n')[0].split(' ')[-1])
                if score > best_score:
                    best_prefix, best_score = prefix, score
        self.log.info('best topology %s %.2f' % (best_prefix, best_score))

        # set to file path
        self._best_tree = f'{self._dir / best_prefix}.raxml.bestTree'

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
                                   self._dir / prefix, self.cpus), mode='bootstrap')
            active_threads.append(thread)
            thread.start()
        self.log.info('active BS Threads: %s' % ' '.join(self._prefixes))

        # finish threads
        for thread in active_threads:
            thread.join()
            self.log.info('finished %s' % thread.prefix)

        self.log.debug('computing supports')
        # concat all bootstrap trees
        all_bs_trees = self._dir / 'all_bs_trees.nwk'
        with open(all_bs_trees, 'w') as fh:
            for prefix in self._prefixes:
                try:
                    with open(f'{self._dir / prefix}.raxml.bootstraps', 'r') as bootstraps:
                        fh.write(bootstraps.read() + '\n')
                except FileNotFoundError:
                    self.log.warning('results from %s missing' % prefix)
        self.log.debug('all bootstrap trees in %s' % all_bs_trees)

        run_cmd = f'{self._binary} --support --tree "{self._best_tree}" ' \
                  f'--bs-trees "{all_bs_trees}" --prefix "{self._dir / "_sup"}" ' \
                  f'--threads {self.cpus} --bs-metric fbp,tbe --redo'
        _run_sp(run_cmd)

        # copy best trees with the right support values to parent directory
        shutil.move(self._dir / '_sup.raxml.supportFBP', self.args.final_tree + '_FBP.nwk')
        shutil.move(self._dir / '_sup.raxml.supportTBE', self.args.final_tree + '_TBE.nwk')
        self.log.info('final trees %s_[FBP,TBE].nwk' % self.args.final_tree)

    def _check_msa(self):
        """
        Checks if the MSA can be understood by raxml-ng,
        parses it to binary and estimates effective number of threads.
        If replacing is off, the original MSA will be used.
        :return: number of parallel instances possible on this machine,
        number of CPUs recommended by raxml-ng for each one, and new path to MSA.
        """
        # prep paths for checking, reducing and parsing MSA
        paths = [self._dir / prefix for prefix in ['_chk', '_red']]

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
            self.log.error(f'error in MSA check, see {paths[0]}.raxml.log')
            exit(1)
        self.log.info('MSA check ok')

        # if any duplicates were found, a reduced alignment was written. ignore if not in replace mode
        msa = f'{paths[0]}.raxml.reduced.phy' if len(keepers) > 0 else self.msa

        # now parse, and use reduced alignment if duplicates found. resulting .rba will be smaller
        arg = f'{self._binary} --msa "{msa}" --parse --model {self.args.evomodel} --prefix "{paths[1]}"'
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
            return runs, cpus, f'{paths[1]}.raxml.rba'
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
            log.exception('ML user interrupt')
        else:
            log.exception('ML process failed, returned %s\n' % str(e.returncode)
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
