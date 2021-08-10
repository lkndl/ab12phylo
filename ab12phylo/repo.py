import re
import sys
from argparse import Namespace
from pathlib import Path
from time import time

BASE_DIR = Path(__file__).resolve().parents[1]
TOOLS = BASE_DIR / 'ab12phylo_cmd' / 'tools'
EXE = '.exe' if sys.executable.endswith('.exe') else ''

PATHS = Namespace(**{
    'raw_msa': Path('Trim') / 'raw_msa.fasta',
    'import_msa': Path('import') / 'import_raw_msa.fasta',
    'msa': 'msa.fasta',
    'msa_anno': 'msa_annotated.fasta',
    'missing': 'missing_samples.tsv',
    'tsv': 'metadata.tsv',
    'xml': Path('BLAST') / 'local_blast+_result.xml',
    'www_xml': Path('BLAST') / 'online_blast_result.xml',
    'missing_fasta': Path('BLAST') / 'missing.fasta',
    'bad_seqs': 'bad_seqs.tsv',
    'preview': Path('Trim') / 'trim_preview.png',
    'cbar': Path('Trim') / 'colorbar.png',
    'gbar': Path('Trim') / 'gene_colorbar.png',
    'left': Path('Trim') / 'msa_gbl_pre.png',
    'right': Path('Trim') / 'msa_gbl_post.png',
    'phylo_msa': Path('Trim') / 'msa_phylo.png',
    'tbe': 'tree_TBE.nwk', 'tben': 'tree_TBE_annotated.nwk',
    'fbp': 'tree_FBP.nwk', 'fbpn': 'tree_FBP_annotated.nwk',
    'icon_path': BASE_DIR / 'ab12phylo' / 'files' / 'favi.png',
    'modified_tree': 'modified_tree.nwk',
    'cmd_config': BASE_DIR / 'ab12phylo_cmd' / 'config' / 'conf.cfg'
})


def toalgo(c, suffix=False):
    algo = {'MAFFT': 'mafft', 'Clustal Omega': 'clustalo', 'MUSCLE': 'muscle', 'T-Coffee': 'tcoffee',
            'RAxML-NG': 'raxml-ng', 'IQ-Tree2': 'iqtree2', 'FastTree': 'FastTree'}[c]
    if suffix and sys.platform in ['win32', 'cygwin']:
        algo += '.exe'
    return algo


toname = lambda c: {'raxml-ng': 'RAxML-NG', 'iqtree2': 'IQ-Tree2',
                    'iqtree2.exe': 'IQ-Tree2'}.get(c, 'unknown')
NUCLEOTIDES = ['A', 'C', 'G', 'T', 'N', 'else', '-', ' ', 'S', '?']
toint_map = dict(zip(NUCLEOTIDES, range(len(NUCLEOTIDES))))
toint = lambda c: toint_map.get(c, 5)
seqtoint = lambda s: list(map(toint, s))
toseq = lambda c: NUCLEOTIDES[c]
inttoseq = lambda s: ''.join(list(map(toseq, s)))
togray_map = dict(zip(NUCLEOTIDES, seqtoint('N') * 5 + list(range(5, 9))))
togray = lambda c: togray_map.get(c, 5)
seqtogray = lambda s: list(map(togray, s))

rgx = lambda c: re.compile(r'|'.join(
    ['.*' + re.escape(word.strip()) + '.*' for word in
     c.strip().strip(',').split(',')]))

tohex = lambda c: '#' + ''.join([('0' + hex(min(255, int(round(a * 256))))[2:])[-2:].upper() for a in c])
torgba = lambda a: tuple(round(int(a[1 + i * 2: 3 + i * 2], 16) / 255, 2) for i in range(len(a) // 2))

technicolor = {
    'A': (0.92, 1, 0.4, 1),
    'C': (0.46, 1, 0.44, 1),
    'G': (0.16, 0.44, 0.8, 1),
    'T': (1, 0.47, 0.66, 1),
    'N': (0.84, 0.84, 0.84, 0.6),
    'else': (0, 0, 0, 1),
    '-': (1, 1, 1, 0),
    ' ': (1, 1, 1, 0),
    'S': (1, 1, 1, 0),
    '?': (0, 0, 0, 0)
}

green_blue = dict(zip(NUCLEOTIDES, [
    (.56, 1, .14, 1),
    (.16, .44, .8, 1),
    (.06, 1, .75, 1),  # rgb(15,255,192)
    (.92, 1, .7, 1),
    (.94, .94, .94, .6),
    (0, 0, 0, 1),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (0, 0, 0, 0)]))

viridis = dict(zip(NUCLEOTIDES, [
    (.48, .11, .44, 1),
    (.22, .55, .77, 1),
    (.1, .9, .51, 1),
    (1, .9, 0, 1),
    (.94, .94, .94, .6),
    (0, 0, 0, 1),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (0, 0, 0, 0)]))

clustalx = dict(zip(NUCLEOTIDES, [
    (.9, .2, .1, 1),  # clustal-red
    (.1, .5, .9, 1),  # clustal-blue
    (.9, .6, .3, 1),  # clustal-orange
    (.1, .8, .1, 1),  # clustal-green
    (.4, .4, .4, 1),  # clustal-dark-gray
    (0, 0, 0, 1),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (0, 0, 0, 0)]))

blue_pink = dict(zip(NUCLEOTIDES, [
    (.4, .65, .96, 1),
    (.14, .05, .8, 1),
    (.99, .33, .93, 1),
    (1, .91, .96, 1),
    (.94, .94, .94, .6),
    (0, 0, 0, 1),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (1, 1, 1, 0),
    (0, 0, 0, 0)]))

colors = list(map(tohex, map(technicolor.get, NUCLEOTIDES)))
USER = 'ab12phylo@gmail.com'
SEP = '??'
ALPHA = .25
DPI = 300
BUF_SIZE = 128 * 1024

DOWNLOAD_TIMEOUT = 3600  # an hour at most for db download

# support value colors
blue = (.2, .5, .75)
red = (1, 0, .3)
dark_red = (.8, 0, .25)
black = (.2, .2, .2)

# seaborn rocket palette
rocket = [blue] + [(0.1237, 0.0717, 0.1822), (0.2452, 0.1049, 0.2639),
                   (0.3809, 0.1206, 0.3250), (0.5172, 0.1179, 0.3545),
                   (0.6582, 0.0955, 0.3536), (0.7965, 0.1050, 0.3106),
                   (0.8949, 0.2178, 0.2531), (0.9429, 0.3754, 0.2636),
                   (0.9592, 0.5330, 0.3748), (0.9644, 0.6702, 0.5150),
                   (0.9689, 0.7980, 0.6851)]

version_regex = re.compile('\\.[\\d]+$')


def inc_priv_timestamp():
    return str(time())[6:17]


help = {0: '<b>Welcome to AB12PHYLO!</b>      To show or hide this help, press '
           '<b>F1</b> or <b>Ctrl+H</b>.\nPlease define your dataset of ABI trace '
           'files or sequence data in <i>.fasta</i> format here. For automatic '
           'mapping of reference data from different genes to the same strain, it is '
           'recommended to use GenBank <i>.fasta</i> records such as <a href="https://'
           'www.ncbi.nlm.nih.gov/nuccore/AF347033.1?report=fasta">AF347033.1</a>. '
           'If wellsplates are used in conjunction with sequence data for multiple '
           'genes, please make sure the plate layouts are identical.',
        1: 'Parse sample ID / well coordinates, gene, orientation and wellsplate ID '
           'from the file name using <a href="https://regex101.com/r/Yulwlf/8" title='
           '"development stage examples">Regular Expressions</a> with capturing groups. '
           'For <i>.fasta</i> files with multiple records, the sequence ID is used. '
           'Fields can be edited manually.',
        2: 'ABI trace files can be trimmed to remove low-quality regions here, based '
           'on phred quality scores. The maximum score is 60, indicating the probability '
           'of an incorrect base call is 1:10^6. Equivalently, 30 indicates 1:1000.',
        3: 'Create per-gene multiple sequence alignments, which are then '
           'concatenated for ML tree inference. AB12PHYLO can construct '
           'MSAs locally, use the <a href="https://www.ebi.ac.uk/Tools/msa/">'
           'EMBL-EBI</a> API or import previously created data.',
        4: 'Trim alignments using <a href="http://molevol.cmima.csic.es/castresana/'
           'Gblocks.html">Gblocks</a>. Especially for mixed datasets, consistent ends '
           'of the MSAs are important to avoid artificial grouping of sequences with '
           'similar sequencing characteristics.',
        5: 'Optionally include species annotation from a local BLAST+ in a '
           '<a href="https://ftp.ncbi.nlm.nih.gov/blast/db/">pre-compiled</a>'
           ' or custom database, or by importing XML results of a manual web '
           '<a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi">BLAST</a> '
           'for one of the <i>_project_/_gene_/_gene_.fasta</i> files. '
           'Searching an online database such as <a href="'
           'https://www.ncbi.nlm.nih.gov/nucleotide/">NCBI nt</a> via the public '
           'BLAST API is also possible, but should not be the main search strategy. '
           'You can continue while BLAST is running.',
        6: 'Infer a phylogenetic tree using <a href="https://github.com/amkozlov/raxml'
           '-ng/">RAxML-NG</a> in three stages: <b>1:</b> Find the maximum-likelihood '
           'tree for the specified evolutionary model from multiple tree searches. <b>'
           '2:</b> Estimate confidence values for the phylogenetic placement of the taxa '
           'in the ML tree by Bootstrapping. <b>3:</b> '
           'Map the result of <b>2</b> to the tree from <b>1</b> as branch support values.',
        7: 'Plot the resulting ML tree and calculate basic diversity / neutrality metrics '
           'in the dataset. Click inside the visualizations to select samples, and right-'
           'click for tree operations such as rooting, dropping nodes, collapsing or '
           'extracting taxa. Species labels can be edited on the BLAST page, and shared '
           'background color indicates identical sequences.'
        }
