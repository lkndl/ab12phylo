#!/usr/bin/python3

"""
Makes result HTML interactive by enabling motif search. Page is otherwise identical.
"""

import cgi
import cgitb
import pickle
import re
import toytree

from os import path, makedirs
from lxml import etree, html
from lxml.html import builder as E
from Bio import SeqIO
from dendropy import DnaCharacterMatrix
from dendropy.calculate import popgenstat
from toytree.utils import ToytreeError

cgitb.enable()
# cgi.test()

# get passed form data
form = cgi.FieldStorage()
call = form['submit'].value
motifs = form['motifs'].value
msa = form['msa_path'].value
pickle_file = form['pickle'].value
try:
    exclude = form['exclude'].value
except KeyError:
    exclude = ''
as_files = True if form['out_type'].value == 'files' else False  # or exclude != '' else False


def _process(leaves, ex_motifs, as_files):
    """
    Raise error if no leaves found.
    If exclusion patterns were specified, remove matching motifs.
    If requested or samples were excluded, return file paths.
    :return:
    """
    print('Content-type: text/html\n\n')
    if len(leaves) == 0:
        raise ValueError('no matches')

    excluded = list()
    if ex_motifs != '':
        # compile regex from patterns
        queries = ['.*' + re.escape(word.strip()) + '.*' for word in ex_motifs.split(',')]
        regex = re.compile(r'|'.join(queries))

        for ex in list(filter(regex.match, leaves)):
            excluded.append(ex)
            leaves.remove(ex)
        # as_files = True if len(excluded) > 0 or as_files else False

    # cut to sample IDs
    leaves = [leaf.split(' ')[0] for leaf in leaves]
    if not as_files:
        return leaves, excluded, None

    # translate to file paths
    with open('metadata.tsv', 'r') as tsv:
        sam_to_file = {leaves.index(sample): _file for sample, _file
                       in [line.split('\t')[1:3] for line in tsv.readlines()] if sample in leaves}
    files = [_file for sample, _file in sorted(sam_to_file.items())]
    return leaves, excluded, files


# read toytree from pickle file
with open(pickle_file, 'rb') as pickle_fh:
    tree = pickle.load(pickle_fh)

# fetch sample IDs
try:
    makedirs('queries', exist_ok=True)

    # compile regex from motifs. empty motifs is caught before
    queries = ['.*' + re.escape(word.strip()) + '.*' for word in motifs.split(',')]
    regex = re.compile(r'|'.join(queries))

    if call == 'subtree':
        mrca = tree.get_mrca_idx_from_tip_labels(regex=regex)
        leaves, excluded, files = _process(tree.get_tip_labels(mrca), exclude, as_files)
        txt = '\n'.join(reversed(files)) if as_files else '\n'.join(reversed(leaves))
        mrca = '%d match%s, subtree MRCA idx: %s' % (len(leaves), 'es' if len(leaves) > 1 else '', mrca)
    elif call == 'match':
        leaves, excluded, files = _process(list(filter(regex.match, tree.get_tip_labels())), exclude, as_files)
        txt = '\n'.join(files) if as_files else '\n'.join(leaves)
        mrca = '%d match%s' % (len(leaves), 'es' if len(leaves) > 1 else '')
    else:
        raise ValueError('illegal call')

except (ToytreeError, ValueError) as ex:
    txt = str(ex)
    mrca = ''
    leaves = []

# save result in file
_path = path.join('queries', call + '_' + motifs + '_' + exclude + '.txt')
with open(_path, 'w') as txt_fh:
    txt_fh.write(txt + '\n')

if len(leaves) > 1:
    # get filtered seqs from MSA
    records = {record.id.split(' ')[0]: record.seq
               for record in SeqIO.parse(msa, 'fasta')
               if record.id.split(' ')[0] in leaves}

    # calc diversity stats
    dna = DnaCharacterMatrix.from_dict(records)
    S = popgenstat.num_segregating_sites(dna)
    pi = popgenstat.nucleotide_diversity(dna)
    theta = popgenstat.wattersons_theta(dna)
    haplo = len(set(records.values())) if S != 0 else 1
    try:
        TD = popgenstat.tajimas_d(dna)
    except ZeroDivisionError:
        TD = float('inf')

    # write raw html table
    table = '<table border="1" class="dataframe">' \
            '<thead><tr style="text-align: justify;"><th>Segregating Sites</th><th>Nucleotide diversity π' \
            '</th><th>Watterson\'s θ</th><th>Tajima\'s D</th><th>Haplotypes</th></tr></thead>' \
            '<tbody><tr><td>%d</td><td>%.5f</td><td>%.5f</td><td>%.5f</td><td>%d</td></tr></tbody></table><hr/>' \
            % (S, pi, theta, TD, haplo)
else:
    table = ''

# read in original page
ET = etree.parse('result.html', parser=etree.HTMLParser())
port = int(ET.findall('//meta[@name="port"]')[0].get('content'))

# fetch div of motif search results
div = ET.findall('//div[@id="motif_result"]')[0]
div.clear()

# insert matches, mrca info, file link and stats into HTML
div.append(E.PRE(E.CODE(txt)))
div.append(E.P(mrca, E.SPAN('see ', E.A('file', href='http://localhost:%s/%s' % (port, _path)), style='float:right')))
div.append(etree.fromstring('<div>%s</div>' % table))

# edit logo and icon path
ET.findall('//img[@id="logo"]')[0].set('src', 'http://localhost:%d/.img/favicon.png' % port)
ET.findall('//link[@id="favicon"]')[0].set('href', 'http://localhost:%d/.img/favicon.ico' % port)

# re-set metadata.tsv link
ET.findall('//a[@id="tsv_link"]')[0].set('a', '../metadata.tsv')
# paste into input fields
ET.findall('//input[@id="motifs"]')[0].set('value', motifs)
ET.findall('//input[@id="exclude"]')[0].set('value', exclude)

# convert to raw text
html = html.tostring(ET, encoding='unicode')

# create a markup for selected nodes
markup = ';font-weight: bold; font-style: italic; text-decoration: underline 0.15rem">%s'
no_markup = '">%s'

# remove old mark-up
for leaf in [tip.split(' ')[0] for tip in tree.get_tip_labels()]:
    html = html.replace(markup % leaf, no_markup % leaf)

# insert new mark-up
for leaf in leaves:
    html = html.replace(no_markup % leaf, markup % leaf)

print('Content-type: text/html\n\n')
print(html)
