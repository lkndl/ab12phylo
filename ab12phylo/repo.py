from argparse import Namespace
from pathlib import Path
import re
from time import time

BASE_DIR = Path(__file__).resolve().parents[1]
TOOLS = Path(__file__).resolve().parents[1] / 'ab12phylo_cmd' / 'tools'

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
    'RAxML': Path('raxml-ng_v1.0.1_linux_x86_64') / 'raxml-ng',
    'tbe': 'tree_TBE.nwk', 'tben': 'tree_TBE_annotated.nwk',
    'fbp': 'tree_FBP.nwk', 'fbpn': 'tree_FBP_annotated.nwk',
    'icon_path': str(Path(__file__).resolve().parent / 'files' / 'favi.png'),
    'modified_tree': 'modified_tree.nwk'
})

algos = {'MAFFT': 'mafft', 'Clustal Omega': 'clustalo', 'MUSCLE': 'muscle', 'T-Coffee': 'tcoffee',
         'RAxML-NG': 'raxml-ng', 'IQ-Tree': 'iqtree2', 'FastTree': 'FastTree'}
toalgo = lambda c: algos[c]
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

tohex = lambda c: '#' + ''.join([(hex(min(255, int(round(a * 256))))[2:] + '0')[:2].upper() for a in c])

KXLIN = {
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

colors = list(map(tohex, map(KXLIN.get, NUCLEOTIDES)))
USER = 'leo.kaindl@tum.de'
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
inc_priv_timestamp = lambda c: str(time()).replace('.', '')[6:]

help = {0: '<b>Welcome to AB12PHYLO!</b> \nPlease define your dataset of ABI trace '
           'files or sequence data in <i>.fasta</i> format here. For automatic '
           'mapping of reference data from different genes to the same strain, it is '
           'recommended to use GenBank <i>.fasta</i> records such as <a href="https://'
           'www.ncbi.nlm.nih.gov/nuccore/AF347033.1?report=fasta">AF347033.1</a>. '
           'If wellsplates are used in conjunction with sequence data for multiple '
           'genes, please make sure the plate layouts are identical.',
        1: 'Parse sample IDs / well coordinates, gene, orientation and '
           'wellsplate ID from the respective file name using <a href='
           '"https://regex101.com/r/Yulwlf/8" title="look at some pre-'
           'publication/development stage examples">Regular Expressions'
           '</a> with capturing groups. Fields may be manually corrected.',
        2: 'ABI trace files can be trimmed to remove low-quality regions here, based '
           'on phred quality scores. The maximum score is 60, indicating the probability '
           'of an incorrect base call is 1:10^6. Analogously, 30 indicates 1:1000.',
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
           '<a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi">BLAST</a>. '
           'Searching an online database such as <a href="'
           'https://www.ncbi.nlm.nih.gov/nucleotide/">NCBI nt</a> via the public '
           'BLAST API is also possible, but should not be the main search strategy.',
        6: 'Infer a phylogenetic tree using <a href="https://github.com/amkozlov/raxml'
           '-ng/">RAxML-NG</a> in three stages: <b>1:</b> Find the maximum-likelihood '
           'tree for the specified evolutionary model from multiple tree searches. <b>'
           '2:</b> Estimate confidence values in this tree by Bootstrapping. <b>3:</b> '
           'Map the result of <b>2</b> to the tree from <b>1</b> as branch support values.',
        7: 'Plot the resulting ML tree and calculate basic diversity / neutrality metrics '
           'in the dataset. Click inside the visualizations to select samples, and right-'
           'click for tree operations such as rooting, dropping nodes, collapsing or '
           'extracting taxa.'
        }
