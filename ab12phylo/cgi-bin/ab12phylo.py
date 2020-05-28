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

# read in original page
with open('result.html', 'r') as html_fh:
    html = html_fh.read()

# read toytree.Toytree.ToyTree from pickle file
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
        # get tips and cut to sample IDs
        leaves = [leaf.split(' ')[0] for leaf in tree.get_tip_labels(mrca)]
        if len(leaves) == 0:
            raise ValueError('no matches')
        txt = '\n'.join(list(reversed(leaves)))
        mrca = '%d match%s, subtree MRCA idx: %s' % (len(leaves), 'es' if len(leaves) > 1 else '', mrca)
    elif call == 'match':
        # get tips and cut to sample IDs
        leaves = [leaf.split(' ')[0] for leaf in filter(regex.match, tree.get_tip_labels())]
        if len(leaves) == 0:
            raise ValueError('no matches')
        txt = '\n'.join(leaves)
        mrca = '%d match%s' % (len(leaves), 'es' if len(leaves) > 1 else '')
    else:
        raise ValueError('illegal call')

except (ToytreeError, ValueError) as ex:
    txt = str(ex)
    mrca = ''
    leaves = []

# save result in file
_path = path.join('queries', call + '_' + motifs + '.txt')
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

    # paste into html table
    table = '<table border="1" class="dataframe">' \
            '<thead><tr style="text-align: justify;"><th>Segregating Sites</th><th>Nucleotide diversity π' \
            '</th><th>Watterson\'s θ</th><th>Tajima\'s D</th><th>Haplotypes</th><tr><thead>' \
            '<tbody><tr><td>%d</td><td>%.5f</td><td>%.5f</td><td>%.5f</td><td>%d</td></tr></tbody></table><hr/>' \
            % (S, pi, theta, TD, haplo)
else:
    table = ''

# edit logo and icon path
html = html.replace('src="cgi-bin/favicon.png"', 'src="http://localhost:8000/.img/favicon.png"')
html = html.replace('href="cgi-bin/favicon.ico"', 'href="http://localhost:8000/.img/favicon.ico"')
html = html.replace('source <a href="metadata.tsv">.tsv</a> contains', 'source <a href="../metadata.tsv">.tsv</a> contains')

# paste into HTML
html = html.replace('<!--insert-->', '<pre id="code_box"><code id="query_result">%s</code></pre>'
                                     '<p>%s<span style="float:right">see '
                                     '<a href="http://localhost:8000/%s">file</a></span></p>%s'
                    % (txt, mrca, _path, table))

# paste into input field
html = html.replace('<input type="text" id="motifs" name="motifs" value=""> &emsp;',
                    '<input type="text" id="motifs" name="motifs" value="%s"> &emsp;' % motifs)

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
