from argparse import Namespace
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[2]
TOOLS = Path(__file__).resolve().parents[2] / 'ab12phylo' / 'tools'

PATHS = Namespace(**{
    'raw_msa': Path('Trim') / 'raw_msa.fasta',
    'import_msa': Path('import') / 'import_raw_msa.fasta',
    'msa': 'msa.fasta',
    'missing': 'missing_samples.tsv',
    'tsv': 'metadata.tsv',
    'xml': Path('BLAST') / 'local_blast+_result.xml',
    'www_xml': Path('BLAST') / 'online_blast_result.xml',
    'missing_fasta': Path('BLAST') / 'missing.fasta',
    'bad_seqs': 'bad_seqs.tsv',
    'preview': Path('Trim') / 'trim_preview.png',
    'cbar': Path('Trim') / 'colorbar.png',
    'left': Path('Trim') / 'msa_gbl_pre.png',
    'right': Path('Trim') / 'msa_gbl_post.png'
})

algos = {'MAFFT': 'mafft', 'Clustal Omega': 'clustalo', 'MUSCLE': 'muscle', 'T-Coffee': 'tcoffee',
         'RAxML-NG': 'raxml-ng', 'IQ-Tree': 'iqtree', 'FastTree': 'FastTree'}
toalgo = lambda c: algos[c]

NUCLEOTIDES = ['A', 'C', 'G', 'T', 'N', 'else', '-', ' ', 'S']
toint_map = dict(zip(NUCLEOTIDES, range(len(NUCLEOTIDES))))
toint = lambda c: toint_map.get(c, 5)
seqtoint = lambda s: list(map(toint, s))
togray_map = dict(zip(NUCLEOTIDES, seqtoint('N') * 5 + list(range(5, 9))))
togray = lambda c: togray_map.get(c, 5)
seqtogray = lambda s: list(map(togray, s))

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
    'S': (1, 1, 1, 0)}

colors = list(map(tohex, map(KXLIN.get, NUCLEOTIDES)))
USER = 'leo.kaindl@tum.de'
SEP = '?*?'
ALPHA = .25
DPI = 300
BUF_SIZE = 128 * 1024
