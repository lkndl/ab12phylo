#!/usr/bin/python3

"""
Makes result HTML interactive by enabling motif search.
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


def _exclude(ex_motifs, leaves):
    # filter leaves with exclusion motifs
    if ex_motifs != '':
        # compile exclusion regex.
        queries = ['.*' + re.escape(word.strip()) + '.*' for word in ex_motifs.split(',')]
        if '.*.*' in queries and len(queries) > 1:
            queries.remove('.*.*')
        regex = re.compile(r'|'.join(queries))

        # drop excluded tips
        for dropped_tip in list(filter(regex.match, leaves)):
            leaves.remove(dropped_tip)
    return leaves


def _to_files(leaves):
    # translate to file paths
    with open('metadata.tsv', 'r') as tsv:
        sam_to_file = {leaves.index(sample): _file for sample, _file
                       in [line.split('\t')[1:3] for line in tsv.readlines()] if sample in leaves}
    return [_file for sample, _file in sorted(sam_to_file.items())]


cgitb.enable()
# cgi.test()
print('Content-type: text/html\n\n')

# get form data
form = cgi.FieldStorage()
call = form['submit'].value
motifs = form['motifs'].value
as_files = True if form['out_type'].value == 'files' else False
try:
    ex_motifs = form['exclude'].value
except KeyError:
    ex_motifs = ''
with open(form['pickle'].value, 'rb') as pickle_fh:
    tree = pickle.load(pickle_fh)

try:
    makedirs('queries', exist_ok=True)

    # compile motifs regex. empty motifs is caught before
    queries = ['.*' + re.escape(word.strip()) + '.*' for word in motifs.split(',')]
    regex = re.compile(r'|'.join(queries))

    # find matching tips
    leaves = list(filter(regex.match, tree.get_tip_labels()))

    # filter with exclusion pattern
    leaves = _exclude(ex_motifs, leaves)

    if len(leaves) == 0:
        raise ValueError('no matches')

    if call == 'subtree':
        mrca = tree.get_mrca_idx_from_tip_labels(names=leaves)
        leaves = tree.get_tip_labels(mrca)
        leaves = _exclude(ex_motifs, leaves)
        leaves = [leaf.split(' ')[0] for leaf in leaves]
        files = _to_files(leaves)
        txt = '\n'.join(reversed(files)) if as_files else ', '.join(reversed(leaves))
        mrca = '%d match%s, subtree MRCA idx: %s' % (len(leaves), 'es' if len(leaves) > 1 else '', mrca)
    elif call == 'match':
        leaves = [leaf.split(' ')[0] for leaf in leaves]
        files = _to_files(leaves)
        txt = '\n'.join(files) if as_files else ', '.join(leaves)
        mrca = '%d match%s' % (len(leaves), 'es' if len(leaves) > 1 else '')
    else:
        raise ValueError('illegal call')

except (ToytreeError, ValueError) as ex:
    txt = str(ex)
    mrca = ''
    leaves = []

# save result in file
_path = path.join('queries', call + '_' + motifs + '_' + ex_motifs + '.txt')
with open(_path, 'w') as txt_fh:
    txt_fh.write(txt.replace(', ', '\n') + '\n')

# diversity stats
if len(leaves) > 1:
    seqs = {record.id.split(' ')[0]: record.seq
            for record in SeqIO.parse(form['msa_path'].value, 'fasta')
            if record.id.split(' ')[0] in leaves}

    dna = DnaCharacterMatrix.from_dict(seqs)
    S = popgenstat.num_segregating_sites(dna)
    pi = popgenstat.nucleotide_diversity(dna)
    theta = popgenstat.wattersons_theta(dna)
    haplo = len(set(seqs.values())) if S != 0 else 1
    try:
        TD = popgenstat.tajimas_d(dna)
    except ZeroDivisionError:
        TD = float('inf')

    # write raw html table
    table = '<table border="1" class="dataframe">' \
            '<thead><tr style="text-align: justify;"><th>Segregating Sites</th><th>Nucleotide diversity π' \
            '</th><th>Watterson\'s θ</th><th>Tajima\'s D</th><th>Haplotypes</th></tr></thead>' \
            '<tbody><tr><td>%d</td><td>%.5f</td><td>%.5f</td><td>%.5f</td><td>%d</td></tr></tbody></table>' \
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
ET.xpath('//img[@id="logo"]')[0].set('src', 'http://localhost:%d/.img/favicon.png' % port)
ET.xpath('//img[@id="msa_viz"]')[0].set('src', 'http://localhost:%d/rectangular_msa.png' % port)
ET.xpath('//link[@id="favicon"]')[0].set('href', 'http://localhost:%d/.img/favicon.ico' % port)

# re-set metadata.tsv link
ET.xpath('//a[@id="tsv_link"]')[0].set('a', '../metadata.tsv')
# paste into input fields
ET.xpath('//input[@id="motifs"]')[0].set('value', motifs)
ET.xpath('//input[@id="exclude"]')[0].set('value', ex_motifs)

# convert to raw text
html = html.tostring(ET, encoding='unicode')

# create a markup for selected nodes
markup = ';font-weight: bold; font-style: italic; text-decoration: underline 0.15rem">%s'
no_markup = '">%s'
no_markup2 = '">%s '
no_markup3 = '">%s<'

# remove old mark-up
for leaf in [tip.split(' ')[0] for tip in tree.get_tip_labels()]:
    html = html.replace(markup % leaf, no_markup % leaf)

# insert new mark-up
for leaf in leaves:
    html = html.replace(no_markup2 % leaf, markup % leaf + ' ')
    html = html.replace(no_markup3 % leaf, markup % leaf + '<')

print(html)
