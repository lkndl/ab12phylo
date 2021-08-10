# 2020 Leo Kaindl

"""
The AB12PHYLO command line interface defines possible options and parses valid arguments
from user input via :class:`sys.argv`, supplemented by the :file:`config/config.yaml` file.
Arguments will be saved as an :class:`argparse.Namespace` object and directly accessed by the
:module:`main` module. Additionally, this module initiates logging.
"""

import argparse
import configparser
import logging
import os
import random
import shutil
import sys
import zipfile
from os import path
from pathlib import Path

import requests
import yaml

from ab12phylo_cmd import phylo
from ab12phylo_cmd.filter import fetch_non_python_tools
from ab12phylo_cmd.__init__ import __version__, __date__


class parser(argparse.ArgumentParser):
    """
    The command line parser of the package.
    Run `ab12phylo-cmd -h` to display all available options.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(prog='ab12phylo-cmd')

        # if empty commandline, show info and help:
        if not args[0]:
            print('AB12PHYLO commandline version %s built on %s\n' % (__version__, __date__))
            args = (['-h'],) + args[1:]

        # [modes]
        mod = parser.add_argument_group(self, title='RUN MODES')
        parts = mod.add_mutually_exclusive_group()
        parts.add_argument('-p1', '--prepare', action='store_true',
                           help='run first part of ab12phylo-cmd, including BLAST but excluding RAxML-NG/IQ-Tree.')
        parts.add_argument('-p2', '--finish', action='store_true',
                           help='run second part of ab12phylo-cmd, beginning with RAxML-NG/IQ-Tree.')
        parts.add_argument('-px', '--add_xml', action='store_true',
                           help='after -p1 run; only read BLAST results. Pass file via -xml.')
        mod.add_argument('-viz', '--visualize', action='store_true',
                         help='invoke ab12phylo-visualize by appending ab12phylo-cmd command.')
        mod.add_argument('-view', '--view', action='store_true',
                         help='invoke ab12phylo-view by appending ab12phylo-cmd command.')

        # [I/O]
        ion = parser.add_argument_group(self, title='FILE I/O')
        ion.add_argument('-dir', '--dir',
                         help='output directory. Defaults to \'./results\'')
        ion.add_argument('-g', '--genes', nargs='+',
                         help='gene(s) to be considered; first argument defines gene for species annotation. '
                              'If set, only ABI traces with a filename matching one of these patterns will be read.')
        ion.add_argument('-abi', '--abi_dir',
                         help='root directory of all ABI trace files. Defaults to current working directory',
                         type=lambda arg: arg if path.isdir(arg) else self.error(
                             '%s: invalid ABI trace file directory' % arg))
        ion.add_argument('-abiset', '--abi_set',
                         help='whitelist file defining subset of ABI traces for the analysis. '
                              'Files must be in or below provided \'--abi_dir\' directory.',
                         type=lambda arg: arg if path.isfile(arg) else self.error(
                             '%s: invalid whitelist file' % arg))
        ion.add_argument('-sampleset', '--sample_set',
                         help='whitelist file defining subset of sample IDs for the analysis. '
                              'Different versions of a sample will be included.',
                         type=lambda arg: arg if path.isfile(arg) else self.error(
                             '%s: invalid sample whitelist' % arg))
        ion.add_argument('-csv', '--csv_dir',
                         help='root directory of .csv files with well-to-isolate coordinates.',
                         type=lambda arg: arg if path.isdir(arg) else 'ignore')
        ion.add_argument('-r2', '--regex_csv', help='RegEx to parse the wellsplate number from a .csv filename. '
                                                    'Use a single capturing group, and double quotes in bash.')
        rxp = ion.add_mutually_exclusive_group()
        rxp.add_argument('-r1', '--regex_abi', help='RegEx to parse plate number, gene name and well position '
                                                    'from an .ab1 filename in that order. Use double quotes in bash.')
        rxp.add_argument('-r3', '--regex_3', nargs=3,
                         help='alternative to --regex_abi. 3 regular expressions to parse plate number, gene name and '
                              'well from an .ab1 filename. Provide the regular expressions in that order without commas'
                              ', but with double quotes from bash. Use a single capturing group in each regex.')
        ion.add_argument('-r4', '--regex_rev', help='RegEx to identify reverse reads from their .ab1 filename.')

        refs = ion.add_mutually_exclusive_group()
        refs.add_argument('-rf', '--ref', nargs='+',
                          help='optional paths of .fasta-formatted files containing reference sequences. '
                               'Files will be matched to genes by order if --genes is set, otherwise by filename.',
                          type=lambda arg: arg if path.isfile(arg) else self.error(
                              'invalid file path(s) with .fasta-formatted reference sequences:\n%s' % arg))
        refs.add_argument('-rd', '--ref_dir', type=self._valid_ref_dir,
                          help='directory of .fasta files with reference sequences. Files will be matched to genes '
                               'by their filename. Set only one option from {--ref, --ref_dir}.')

        # [quality]
        qal = parser.add_argument_group(self, title='QUALITY')
        qal.add_argument('-qal', '--min_phred', type=int,
                         help='minimal phred quality score to define \'good\' bases in ABI trace files.')
        qal.add_argument('-bad', '--bad_stretch', type=int,
                         help='number of consecutive \'bad bases\': Any sequence of bases in an ABI trace file '
                              'with a phred quality score below the minimum and at least as long as the number '
                              'supplied here will be replaced by a sequence of Ns of equal length.')
        qal.add_argument('-end', '--end_ratio', type=self._valid_end_ratio,
                         help='defines a \'good end\' of a sequence in an ABI trace file for trimming. '
                              'Enter as "<int>/<int>".')

        # [BLAST]
        bla = parser.add_argument_group(self, title='BLAST')
        skips = bla.add_mutually_exclusive_group()
        skips.add_argument('-local', '--no_remote', action='store_true',
                           help='NCBI BLAST API queries are de-prioritized very quickly. Set this flag to skip online '
                                'nucleotide BLAST for seqs missing from the local database.')
        skips.add_argument('-none', '--no_BLAST', action='store_true',
                           help='skip BLAST entirely.')
        skips.add_argument('-remote', '--no_local', action='store_true',
                           help='BLAST only in remote database; discouraged.')
        skips.add_argument('-xml', '--BLAST_xml', nargs='+', help='supply available BLAST results as .XML files.')
        bla.add_argument('-db', '--db', help='BLAST+ database to use. Will attempt to download or update this db from '
                                             'https://ftp.ncbi.nlm.nih.gov/blast/db/ unless provided via -dbpath.')
        bla.add_argument('-dbpath', '--dbpath',
                         help='path to directory with a BLAST+ database. Set this option for a user-created database '
                              'or if ab12phylo-cmd is not allowed FTP access. You might have to define -db as well.',
                         type=lambda arg: arg if path.isdir(arg) else self.error(
                             'invalid path to directory containing BLAST database'))
        bla.add_argument('-remotedb', '--remote_db',
                         help='NCBI database to search for sequences not found locally. Default is \'nt\' for DNA.')

        # [MSA]
        msa = parser.add_argument_group(self, title='MSA')
        msa.add_argument('-algo', '--msa_algo', choices=['clustalo', 'mafft', 'muscle', 't_coffee'],
                         help='select an algorithm to build the Multiple Sequence Alignment. Default is MAFFT.')
        msa.add_argument('-gbl', '--gblocks', choices=['skip', 'relaxed', 'balanced', 'default', 'strict'],
                         help='activate/set MSA trimming with Gblocks.')

        phy = parser.add_argument_group(self, title='ML TREE INFERENCE')
        phy.add_argument('-tool', '--ml_tool', choices=['raxml-ng', 'iqtree2'],
                         help='select a tool to re-construct a phylogenetic tree. RAxML-NG is the default '
                              'for Linux, but not available for Windows.')
        phy.add_argument('-st', '--start_trees', type=self._valid_start_trees,
                         help='numbers of starting trees for raxml-ng / iqtree2 tree inference: '
                              '[<int random trees>,<int parsimony-based trees>].')
        phy.add_argument('-bst', '--bootstrap', type=self._valid_bootstrap,
                         help='number of bootstrap trees. The MRE bootstrap convergence test in RAxML-NG may stop '
                              'iteration before the passed number is reached, but this is unlikely. See '
                              'https://doi.org/10.1089/cmb.2009.0179')
        phy.add_argument('-uf', '--ultrafast', action='store_true',
                         help='for tree inference with IQ-Tree 2, ultrafast bootstrapping can be enabled. '
                              'If a value lower than 1000 is passed to `--bootstrap`, the number of iterations will '
                              'automatically be set to this minimum accepted value.')
        mod = phy.add_mutually_exclusive_group()
        mod.add_argument('-evomodel', '--evomodel', help='evolutionary model for tree inference. The default '
                                                         'is GTR+I+G, and no checks are performed!')
        mod.add_argument('-findmodel', '--findmodel', action='store_true',
                         help='FOR IQ-TREE2 ONLY: infer an evolutionary model using ModelFinder before '
                              're-constructing a tree. https://doi.org/10.1038/nmeth.4285')
        phy.add_argument('-s', '--seed', type=int,
                         help='seed value for reproducible tree inference results. Will be random if not set.')
        phy.add_argument('-metric', '--metric', choices=['TBE', 'FBP'],
                         help='bootstrap support metric: Either Felsenstein Bootstrap Proportions (FBP) '
                              'or Transfer Bootstrap Expectation (TBE). https://doi.org/10.1038/s41586-018-0043-0')

        # [visualize]
        viz = parser.add_argument_group(self, title='VISUALIZATION')
        viz.add_argument('result_dir', nargs='?',
                         help='ONLY FOR AB12PHYLO-VISUALIZE/-VIEW: path to earlier run. Pass without keyword!')
        viz.add_argument('-msa_viz', '--msa_viz', nargs='*', choices=['pdf', 'png'],
                         help='also render a rectangular tree with MSA color matrix in chosen format(s). '
                              'Takes some extra time.')
        viz.add_argument('-threshold', '--threshold', type=float,
                         help='limit between 0 and 1 for support value color switch.')
        viz.add_argument('-out_fmt', '--out_fmt', nargs='+', choices=['pdf', 'png', 'svg'],
                         help='output file format for tree graphics.')
        viz.add_argument('-md', '--min_dist', type=float,
                         help='minimal phylogenetic distance between any pair of samples for RAxML-NG.')
        viz.add_argument('-mpd', '--min_plot_dist', type=float,
                         help='minimal artificial distance of a node in the visualization to its parent.')
        viz.add_argument('-drop', '--drop_nodes', nargs='+', type=int,
                         help='drop the node(s) with this idx and its descendants from the tree. '
                              'Use the motif subtree search to find the node or MRCA idx of the target.')
        viz.add_argument('-replace', '--replace_nodes', nargs='+', type=int,
                         help='replace the node(s) with this idx and its descendants with a placeholder.')
        viz.add_argument('-root', '--root', type=int, help='root the tree at the node with this index.')
        viz.add_argument('-supp', '--print_supports', action='store_true',
                         help='print the support values in percent of the optimal value in the rectangular tree.')

        # [popgen]
        gen = parser.add_argument_group(self, title='POPULATION GENETICS')
        gen.add_argument('-gap', '--gap_share', type=self._valid_threshold,
                         help='maximum share of gaps at an MSA site that will be ignored.')
        gen.add_argument('-unk', '--unknown_share', type=self._valid_threshold,
                         help='maximum acceptable proportion of unknown characters at an MSA site.')
        viz.add_argument('-poly', '--poly_allelic', action='store_true',
                         help='accept segregating sites with more than one mutant nucleotide.')

        # [misc]
        level = self.add_mutually_exclusive_group()
        level.add_argument('-i', '--info', action='store_true', help='show some more information in console output.')
        level.add_argument('-v', '--verbose', action='store_true', help='show all runtime information in console.')
        self.add_argument('-c', '--config',
                          default=path.abspath(path.join(path.dirname(__file__), 'config', 'config.yaml')),
                          type=lambda arg: arg if path.isfile(arg) else self.error('%s: invalid .config path'),
                          help='path to .yaml config file with defaults; command line arguments will override them.')
        self.add_argument('-nt', '--max_threads', type=int,
                          help='Limit the number of CPUs to use for AB12PHYLO.')
        self.add_argument('-version', '--version', action='store_true', help='print version information and exit.')
        self.add_argument('-test', '--test', action='store_true', help='Test run.')
        self.add_argument('-q', '--headless', action='store_true',
                          help='do not start a CGI server nor display in browser. For remote use.')
        self.add_argument('-init', '--initialize', action='store_true',
                          help='re-initialize ab12phylo-cmd: Search for existing BLAST+, '
                               'RAxML-NG and IQ-Tree installations, or re-run these.')

        self._initialize(auto=True)

        self.args = self.parse_args(args[0])

        if self.args.version is True:
            sys.exit('ab12phylo: %s' % __version__)

        cfg = Path(__file__).resolve().parent / 'config' / 'conf.cfg'
        if self.args.initialize is True:
            if cfg.is_file():
                cfg.unlink()
            self._initialize(auto=False)
            exit(0)

        # test: switch config + set verbose
        if self.args.test is True:
            print('--TEST RUN--', file=sys.stderr)
            self.args.config = path.abspath(path.join(path.dirname(__file__), 'config', 'test_config.yaml'))
            self.args.verbose = True

        # load additional info from config
        assert self.args.config is not None
        config = yaml.safe_load(open(self.args.config, 'r'))

        # if refs were set manually; and genes as well, or there is only 1 reference -> match refs to genes by order
        by_order = True if self.args.ref is not None \
                           and (self.args.genes is not None or len(self.args.ref) == 1) else False

        config_only = dict()
        # provide defaults for unset options
        for key, val in config.items():
            if key not in self.args:
                # for items without CLI equivalent (avoid bloated CLI)
                config_only[key] = val
            elif self.args.__dict__.get(key) in [None, False]:
                # access the namespace itself without var names -> access dict
                if key in ['abi_dir', 'csv_dir', 'blastdb', 'abi_set', 'sample_set', 'dir'] and val[0] == '$':
                    # deal with relative paths in config for test case
                    val = path.join(path.dirname(path.dirname(__file__)), val[1:])
                elif key == 'ref':
                    # split into list
                    val = [ref.strip() for ref in val.split(',')]
                    # make absolute paths
                    val = [path.join(path.dirname(path.dirname(__file__)), ref[1:])
                           if ref[0] == '$' else ref for ref in val]
                elif key == 'regex_abi' and self.args.regex_3 or key == 'regex_3' and self.args.regex_abi:
                    # do not interfere with user-defined exclusive group
                    continue
                elif key == 'findmodel' and self.args.evomodel or key == 'evomodel' and self.args.findmodel:
                    # as above
                    continue

                self.args.__dict__[key] = val

        # ab12phylo-cmd with --visualize or --view: guess real results path and skip re-parsing
        if len(kwargs) > 0 or self.args.visualize or self.args.view:
            if self.args.visualize:
                print('starting -viz re-plotting run', file=sys.stderr)
            elif self.args.view:
                print('re-view results', file=sys.stderr)
            # look in current working directory and ./results
            found = False
            for outer in [self.args.result_dir, self.args.dir, os.getcwd()]:
                if found:
                    break
                if outer is None:
                    continue
                elif outer == '.':
                    outer = os.getcwd()
                for inner in ['', 'results']:
                    if path.isfile(path.join(outer, inner, 'tree_TBE.nwk')):
                        # filename is hardcoded manually from below.
                        self.args.dir = path.join(outer, inner)
                        found = True
                        break
            if not found:
                sys.exit('Result files not found')

        else:  # normal case

            # move ref file list from ref_dir to ref
            if self.args.ref_dir:
                if type(self.args.ref_dir) == list:
                    self.args.ref = self.args.ref_dir
                else:
                    self.args.ref = self._valid_ref_dir(self.args.ref_dir)
                del self.args.ref_dir

            # check for duplicates in genes and references
            if self.args.genes is not None:
                if len(set(self.args.genes)) < len(self.args.genes):
                    self.error('duplicates in supplied genes')
            if self.args.ref is not None and len(set(self.args.ref)) < len(list(self.args.ref)):
                self.error('duplicates in supplied references')

            # now rebuild a command line and parse it again
            commandline = list()
            for key, val in self.args.__dict__.items():
                if key in ['genes', 'ref', 'regex_3', 'BLAST_xml', 'out_fmt', 'msa_viz'] and val is not None:
                    commandline.append('--%s' % key)
                    [commandline.append(v) for v in val]
                    if key == 'msa_viz' and not val:
                        commandline.append('png')
                elif val not in [None, False, True] or key == 'max_threads':
                    commandline += ['--%s' % key, str(val)]
                elif val is True:
                    commandline.append('--%s' % key)

            self.args = self.parse_args(commandline)

            if sys.platform in ['win32', 'cygwin'] and self.args.ml_tool == 'raxml-ng' and not \
                    (self.args.visualize or self.args.view or self.args.prepare or self.args.add_xml):
                print(f'\033[91mWARNING: RAxML-NG was selected for ML inference on a Windows '
                      f'or cygwin system. This will likely fail.\033[0m Select IQ-Tree 2 by '
                      f'passing \'--ml_tool iqtree2\' or editing the equivalent line in '
                      f'{self.args.config}')
                if not self.args.headless:
                    answer = ''
                    while answer not in {'y', 'yes', 'n', 'no'}:
                        answer = input(f'Continue anyway? [y/n]').lower().strip()
                    if answer in {'n', 'no'}:
                        exit(1)

            # remember type of original ref option
            self.args.by_order = by_order

            # create output directory already
            if self.args.dir is not None:
                os.makedirs(self.args.dir, exist_ok=True)
            else:
                # write results to current directory.
                self.args.dir = ''

            # define a random seed if None was given
            if self.args.seed is None:
                self.args.seed = random.randint(0, 1000)

        # fetch the paths from the config
        cfg_parser = configparser.ConfigParser()
        cfg_parser.read(cfg)
        self.args.cfg = dict()
        if 'Paths' in cfg_parser:
            self.args.cfg.update(dict(cfg_parser['Paths']))

        # set some default values where options would be useless
        self.args.xml = path.join(self.args.dir, 'BLAST', 'local_blast+_result.xml')
        self.args.www_xml = path.join(self.args.dir, 'BLAST', 'online_blast_result.xml')
        self.args.bad_seqs = path.join(self.args.dir, 'bad_seqs.tsv')
        self.args.missing_samples = path.join(self.args.dir, 'missing_samples.tsv')
        self.args.tsv = path.join(self.args.dir, 'metadata.tsv')
        self.args.msa = path.join(self.args.dir, 'msa.fasta')
        self.args.new_msa = path.join(self.args.dir, 'msa_annotated.fasta')
        self.args.mview_msa = path.join(self.args.dir, 'msa_mview.html')
        self.args.topo = path.join(self.args.dir, 'topology_preview.png')
        self.args.missing_fasta = path.join(self.args.dir, 'missing.fasta')
        self.args.final_tree = path.join(self.args.dir, 'tree')
        self.args.annotated_tree = path.join(self.args.dir, 'tree_%s_annotated.nwk' % self.args.metric)
        if self.args.prepare:
            self.args.log = path.join(self.args.dir, 'ab12phylo-p1.log')
        elif self.args.finish:
            self.args.log = path.join(self.args.dir, 'ab12phylo-p2.log')
        elif self.args.add_xml:
            self.args.log = path.join(self.args.dir, 'ab12phylo-px.log')
        else:
            self.args.log = path.join(self.args.dir, 'ab12phylo.log')
        self.args.sep = 'SSSSSSSSSS'  # to visually separate genes in the concat MSA
        # now also load config-only defaults
        self.args.__dict__.update(config_only)

        # switching to visualize and view
        if self.args.visualize:
            self._init_log(self.args.log[:-4] + '-viz.log')
            log = logging.getLogger(__name__)
            log.debug('--AB12PHYLO-VISUALIZE--')
            log.debug(' '.join(args[0]))
            # copy config
            shutil.copy(src=self.args.config, dst=path.join(self.args.dir, 'used_config.yaml'))
            phylo.tree_build(self.args)
            sys.exit(0)

        if self.args.view:
            self.args.headless = False
            self._init_log(self.args.log[:-4] + '-view.log')
            log = logging.getLogger(__name__)
            log.debug('--AB12PHYLO-VIEW--')
            log.debug(' '.join(args[0]))
            # copy config
            shutil.copy(src=self.args.config, dst=path.join(self.args.dir, 'used_config.yaml'))
            phylo.tree_view(self.args.dir)
            sys.exit(0)

        # configure logging:
        if len(kwargs) > 0:
            self._init_log(self.args.log[:-4] + '-view-viz.log')
            log = logging.getLogger(__name__)
            log.debug(' '.join(args[0]))
        else:
            self._init_log(self.args.log)
            log = logging.getLogger(__name__)
            log.debug('--ARGS-- %s' % ' '.join(args[0]))
            log.debug('running AB12PHYLO v%s' % __version__)
            log.info('seed for this run: %s' % self.args.seed)
            if by_order is True:
                log.info('will match references to genes by order')

        # copy config
        shutil.copy(src=self.args.config, dst=path.join(self.args.dir, 'used_config.yaml'))

    def _initialize(self, auto=False):
        """
        Check if ab12phylo-cmd has been run before, i.e. whether the :file:`conf.cfg`
        file is present. If not, search for BLAST+, RAxML-NG and iqtree2 installations
        via shutil. If a tool is not found on the $PATH or outdated, run prompts and
        installations. If the file is not present, always try downloading test data.
        :return:
        """
        _dir = Path(__file__).resolve().parent
        cfg = _dir / 'config' / 'conf.cfg'
        other_cfg = _dir.parent / 'ab12phylo' / 'conf.cfg'
        if cfg.is_file():
            return

        # init temporary console logger
        log = logging.getLogger(__name__)
        log.setLevel(logging.DEBUG)
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.DEBUG)
        sh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        log.addHandler(sh)

        log.info('Initializing ab12phylo-cmd:')

        log.info('Downloading test data ...')
        r = requests.get('https://github.com/lkndl/ab12phylo/wiki/cmd_test_data.zip',
                         stream=True, timeout=12)
        zf = _dir / 'test_data.zip'
        with open(zf, 'wb') as file:
            for chunk in r.iter_content(chunk_size=128):
                file.write(chunk)
        with zipfile.ZipFile(zf, 'r') as zo:
            zo.extractall(zf.with_suffix(''))

        fetch_non_python_tools('-cmd', cfg, other_cfg, _dir / 'tools', log)
        sh.close()
        log.removeHandler(sh)
        if auto:
            print(f'\033[94mInitialized AB12PHYLO, now exiting.\033[0m')
            exit(0)

    def _valid_ref_dir(self, ref_dir):
        """
        Looks for .fasta files in --ref_dir arg, return as list. Trick :later move to --ref flag.
        :return: .fasta reference files as list
        """
        if not path.isdir(ref_dir):
            raise self.error('invalid references directory: %s' % ref_dir)

        ref_files = list()
        for root, dirs, files in os.walk(ref_dir):
            ref_files += [path.join(root, file) for file in files if file.endswith('.fasta')]

        if len(ref_files) == 0:
            raise self.error('no .fasta references found in directory: %s' % ref_dir)
        return ref_files

    def _valid_end_ratio(self, end_ratio):
        """Checks if --end_ratio argument for trimming is in right format and meaningful"""
        try:
            ratio = [int(d) for d in end_ratio.strip().split('/')]
            if len(ratio) == 2 and ratio[0] <= ratio[1]:
                return ratio
            else:
                raise ValueError
        except ValueError:
            raise self.error('invalid end ratio defined: %s' % end_ratio)

    def _valid_start_trees(self, start_trees):
        """Checks if --start_trees argument is in the right format: [int,int]"""
        try:
            start = [int(d) for d in start_trees[1:-1].split(',')]
            if len(start) == 2:
                return start
            else:
                raise ValueError
        except ValueError:
            raise self.error('invalid start trees: %s' % start_trees)

    def _valid_bootstrap(self, bootstrap):
        """Checks if --bootstrap argument is a number > 1"""
        try:
            bootstrap = int(bootstrap)
            if bootstrap > 1:
                return bootstrap
            else:
                raise ValueError
        except ValueError:
            raise self.error('Number of bootstrap trees must be an int > 1')

    def _valid_threshold(self, t):
        """Checks if popgen thresholds are floats in [0:0.5]"""
        try:
            thresh = float(t)
            if thresh < 0 or thresh > 1:
                raise ValueError
            else:
                return thresh
        except ValueError:
            raise self.error('invalid numerical float threshold: %s' % t)

    def _init_log(self, filename):
        """Initializes logging"""
        log = logging.getLogger()
        log.setLevel(logging.DEBUG)

        # init verbose logging to file
        fh = logging.FileHandler(filename=filename, mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s\t%(name)s\t%(message)s',
                                          datefmt='%Y-%m-%d %H:%M:%S'))
        log.addHandler(fh)

        # init shortened console logging
        sh = logging.StreamHandler(sys.stdout)
        if self.args.verbose:
            sh.setLevel(logging.DEBUG)
        elif self.args.info:
            sh.setLevel(logging.INFO)
        else:
            sh.setLevel(logging.WARNING)
        sh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        log.addHandler(sh)
