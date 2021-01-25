from argparse import Namespace
from pathlib import Path
import re

BASE_DIR = Path(__file__).resolve().parents[1]
TOOLS = Path(__file__).resolve().parents[1] / 'ab12phylo' / 'tools'

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
