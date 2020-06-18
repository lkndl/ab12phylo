#!/usr/bin/python3

"""
Makes result HTML interactive by enabling motif search.
"""

import cgi
import cgitb
import itertools

import numpy as np
import os
import pickle
import re
import toytree

from os import path, makedirs
from lxml import etree, html
from lxml.html import builder as E
from Bio import SeqIO
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


def _h(a):
    # harmonic number
    b = 0
    for c in range(1, a):
        b += 1 / c
    return b


def _qh(a):
    # quadratic harmonic number
    b = 0
    for c in range(1, a):
        b += 1 / np.power(c, 2)
    return b


def _diversity_stats(gene_now, records, limits, poly):
    # also accept polyallelic sites as segregating sites

    # translate dictionary of Seqs to numpy int array
    seqs = np.array([seq for seq in records.values()])
    n = len(records)
    bases = ['A', 'C', 'G', 'T', '-']
    weird_chars = np.setdiff1d(np.unique(seqs), bases)
    codes = dict(zip(bases + list(weird_chars), range(len(bases) + len(weird_chars))))
    coded_seqs = np.vectorize(codes.get)(seqs)

    conserved, biallelic, polyallelic, unknown, gaps = 0, 0, 0, 0, 0
    drop = set()

    # iterate over MSA columns, search for segregating sites
    for j in range(coded_seqs.shape[1]):
        col = coded_seqs[:, j]
        if len(set(col)) == 1:
            # singleton
            it = col[0]
            if it > 4:
                unknown += 1
                drop.add(j)
            elif it == 4:
                gaps += 1
                drop.add(j)
            else:
                conserved += 1
        else:  # more than one item
            if max(col) > 4 and len(col[col > 4]) > limits[1] * n:
                unknown += 1
                drop.add(j)
            elif max(col) == 4 and len(col[col >= 4]) > limits[0] * n:
                gaps += 1
                drop.add(j)
            else:
                # no gap, no unknown, more than one
                if len(set(col[col < 4])) > 2:
                    polyallelic += 1
                    if not poly:
                        drop.add(j)
                else:
                    biallelic += 1
    # print('allsites:%d\tconserved:%d\tbiallelic:%d\tpolyallelic:%d\tunknown:%d\tgaps:%d\tpoly:%s\ndrop:%d'
    #       % (coded_seqs.shape[1], conserved, biallelic, polyallelic, unknown, gaps, poly, len(drop)))

    seg_sites = biallelic
    n_sites = coded_seqs.shape[1] - unknown
    if poly:
        seg_sites += polyallelic
    else:
        n_sites -= polyallelic

    if seg_sites == 0:
        return '<tr><td>%s</td><td>0</td><td>--</td><td>--</td><td>--</td><td>--</td>' \
               '<td>--</td><td>%d</td><td>%d</td></tr>' % (gene_now, gaps, unknown)

    # crop to allowed sites
    coded_seqs = coded_seqs[:, list(set(range(coded_seqs.shape[1])) - drop)]

    pi = 0
    for combi in itertools.combinations(range(n), 2):
        pi += np.sum(coded_seqs[combi[0]] != coded_seqs[combi[1]])
    pi = pi * 2 / n / (n - 1)
    pi_per_site = pi / n_sites

    theta_w = seg_sites / _h(n)

    # Tajima's D
    try:
        e1 = 1 / _h(n) * ((n + 1) / (3 * n - 3) - 1 / _h(n))
        b2 = (2 * (n * n + n + 3)) / (9 * n * (n - 1))
        c2 = b2 - ((n + 2) / (n * _h(n))) + (_qh(n) / (_h(n) * _h(n)))
        e2 = c2 / (_h(n) * _h(n) + _qh(n))
        td = (pi - theta_w) / np.sqrt(e1 * seg_sites + e2 * seg_sites * (seg_sites - 1))
    except (ZeroDivisionError, RuntimeWarning) as z:
        td = float('inf')

    haplo = len(np.unique(coded_seqs, axis=0))

    return '<tr><td>%s</td><td>%d</td><td>%.5f</td><td>%.5f</td>' \
           '<td>%.5f</td><td>%.5f</td><td>%d</td><td>%d</td><td>%d</td></tr>' \
           % (gene_now, seg_sites, pi, pi_per_site, theta_w, td, haplo, gaps, unknown)


cgitb.enable()
# cgi.test()
print('Content-type: text/html\n\n')

# get form data
form = cgi.FieldStorage()
call = form['submit'].value
motifs = form['motifs'].value
thresholds = float(form['gap'].value), float(form['unknown'].value)
poly_allowed = True if form['poly'].value == 'True' else False
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
    # per-gene metrics and Tajima's D in sliding-window
    pop_seqs = {record.id.split(' ')[0]: record.seq
                for record in SeqIO.parse(form['msa_path'].value, 'fasta')
                if record.id.split(' ')[0] in leaves}

    # get gene lengths from form
    g_lens = [entry.split(':') for entry in form['g_lens'].value.split('_')]
    pos = 0
    rows = ''

    # calculate diversity/neutrality statistics per gene
    for entry in g_lens:
        gene, glen = entry[0], int(entry[1])
        cut_seqs = {seq_id: seq[pos:pos + glen] for seq_id, seq in pop_seqs.items()}
        rows += _diversity_stats(gene, cut_seqs, thresholds, poly_allowed)
        pos += glen

    # overall if sensible
    if len(g_lens) > 1:
        rows += _diversity_stats('overall', pop_seqs, thresholds, poly_allowed)

    # write raw html table
    table = '<hr/><h4>Neutrality + Diversity Statistics</h4><p>Exclusion threshold for gap sites ' \
            'was <code>%.2f</code>, and <code>%.2f</code> for unknown sites.</p>' \
            '<table border="1" class="dataframe"><thead><tr style="text-align: justify;"><th></th>' \
            '<th># Segregating Sites</th><th>Nucl. diversity π</th><th>π per site</th><th>Watterson\'s θ</th>' \
            '<th>Tajima\'s D</th><th>Haplotypes</th><th>gaps</th><th>unknown</th></tr></thead>' \
            '<tbody>%s</tbody></table>' % (thresholds[0], thresholds[1], rows)
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
