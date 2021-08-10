#!/usr/bin/env python3
# 2021 Leo Kaindl

"""
Makes result HTML interactive by enabling motif search.
Compute neutrality/diversity stats for selected population.
"""

import cgi
import cgitb
import itertools
import sys

import numpy as np
import pickle
import re
import pandas as pd
import toytree

from os import path, makedirs
from lxml import etree, html
from lxml.html import builder as E
from Bio import SeqIO
from toytree.utils import ToytreeError

blue = (0.1960, 0.5333, 0.7411)


def _exclude(ex_motifs, leaves):
    # filter leaves with exclusion motifs
    if ex_motifs != '':
        # compile exclusion regex.
        queries = ['.*' + re.escape(word.strip()) + '.*' for word in ex_motifs.split(',')]
        if '.*.*' in queries and len(queries) > 1:
            queries.remove('.*.*')
        regex = re.compile(r'|'.join(queries))

        # drop excluded tips
        for dropped_tip in list(filter(regex.match, leaves.keys())):
            leaves.pop(dropped_tip)
    return leaves


def _to_files(leaves):
    # translate to file paths
    with open('metadata.tsv', 'r') as tsv:
        # prep an empty list for each leaf
        sam_to_file = {i: list() for i in range(len(leaves))}
        # iterate over the table and for matching IDs append file paths
        [sam_to_file[leaves.index(sample)].append(_file) for sample, _file
         in [line.split('\t')[0:3:2]  # get first and third column
             for line in tsv.readlines()] if sample in leaves]
    return ['\n'.join(files) for sample, files in sorted(sam_to_file.items())]


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
    # translate dictionary of sequences to numpy int array
    seqs = np.array([seq for seq in records.values()])
    n_sequences = len(records)
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
            # only one char at site
            it = col[0]
            if it > 4:
                # all-unknown-site. can happen as we selected a subset of sequences
                unknown += 1
                drop.add(j)
            elif it == 4:
                # all-gap-site
                gaps += 1
                drop.add(j)
            else:
                conserved += 1
        else:  # more than one character in the whole column
            if max(col) > 4 and len(col[col > 4]) > limits[1] * n_sequences:
                # there are too many unknown chars in the column
                unknown += 1
                drop.add(j)
            elif max(col) == 4 and len(col[col >= 4]) > limits[0] * n_sequences:
                # there are too many gaps at the site
                gaps += 1
                drop.add(j)
            else:
                # if gaps or unknown characters are at the site, replace them with the most common nucleotide
                try:
                    col[col >= 4] = np.bincount(col[col < 4]).argmax()
                except ValueError:
                    # there are no nucleotides
                    drop.add(j)
                    continue

                # ordinary site treatment
                num_diff_chars = len(set(col))
                if num_diff_chars > 2:
                    polyallelic += 1
                    if not poly:
                        drop.add(j)
                elif num_diff_chars == 2:
                    biallelic += 1
                else:
                    conserved += 1
    # print('GENE:%s\tallsites:%d\tconserved:%d\tbiallelic:%d\tpolyallelic:%d\tunknown:%d\tgaps:%d\tpoly:%s\ndrop:%d\n'
    #       % (gene_now, coded_seqs.shape[1], conserved, biallelic, polyallelic, unknown, gaps, poly, len(drop)))

    seg_sites = biallelic
    n_sites = coded_seqs.shape[1] - unknown - gaps
    if poly:
        seg_sites += polyallelic
    else:
        n_sites -= polyallelic

    if seg_sites == 0:
        return '<tr><td>%s</td><td>%d</td><td>0</td><td>--</td><td>--</td><td>--</td><td>--</td>' \
               '<td>--</td><td>%d</td><td>%d</td></tr>' % (gene_now, n_sites, gaps, unknown)

    # crop to allowed sites
    coded_seqs = coded_seqs[:, list(set(range(coded_seqs.shape[1])) - drop)]

    k = 0
    # count pairwise differences
    for combi in itertools.combinations(range(n_sequences), 2):
        k += np.sum(coded_seqs[combi[0]] != coded_seqs[combi[1]])

    # binomial coefficient, k = k^ in formula (10) from Tajima1989
    k = k * 2 / n_sequences / (n_sequences - 1)
    pi = k / n_sites

    a1 = _h(n_sequences)
    theta_w = seg_sites / a1 / n_sites  # per site from Yang2014

    # Tajima's D Formula:
    # d = k - seg_sites/a1 = k - theta_W
    # D = d / sqrt(Var(d))
    try:
        a2 = _qh(n_sequences)
        e1 = 1 / a1 * ((n_sequences + 1) / (3 * n_sequences - 3) - 1 / a1)
        b2 = (2 * (n_sequences * n_sequences + n_sequences + 3)) / (9 * n_sequences * (n_sequences - 1))
        c2 = b2 - ((n_sequences + 2) / (n_sequences * a1)) + (a2 / (a1 * a1))
        e2 = c2 / (a1 * a1 + a2)
        tajima = (k - theta_w) / np.sqrt(e1 * seg_sites + e2 * seg_sites * (seg_sites - 1))
    except (ZeroDivisionError, RuntimeWarning) as z:
        tajima = float('inf')

    n_genotypes = len(np.unique(coded_seqs, axis=0))

    return '<tr><td>%s</td><td>%d</td><td>%d</td><td>%.5f</td><td>%.5f</td>' \
           '<td>%.5f</td><td>%.5f</td><td>%d</td><td>%d</td><td>%d</td></tr>' \
           % (gene_now, n_sites, seg_sites, k, pi, theta_w, tajima, n_genotypes, gaps, unknown)


cgitb.enable()
# cgi.test()
print('Content-type: text/html')
print('')

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
    tree_nodes = tree.get_feature_dict(key_attr='name')

    # compile motifs regex. empty motifs is caught before
    queries = ['.*' + re.escape(word.strip()) + '.*' for word in motifs.split(',')]
    regex = re.compile(r'|'.join(queries))

    # find matching tips
    leaf_names = list(filter(regex.match, tree.get_tip_labels()))
    # check if they are references
    leaves = {leaf: tree_nodes[leaf].type == blue for leaf in leaf_names}

    # filter with exclusion pattern
    leaves = _exclude(ex_motifs, leaves)

    if len(leaves) == 0:
        raise ValueError('no matches')

    if call == 'subtree':
        mrca = tree.get_mrca_idx_from_tip_labels(names=leaf_names)
        leaf_names = tree.get_tip_labels(mrca)
        leaves = {leaf: tree_nodes[leaf].type == blue for leaf in leaf_names}
        leaves = _exclude(ex_motifs, leaves)
        leaf_names = [leaf_name.split(' ')[0] if not is_ref else leaf_name
                      for leaf_name, is_ref in leaves.items()]
        files = _to_files(leaf_names)
        txt = '\n'.join(reversed(files)) if as_files else ', '.join(reversed(leaf_names))
        mrca = '%d match%s, subtree MRCA idx: %s' % (len(leaves), 'es' if len(leaves) > 1 else '', mrca)
    elif call == 'match':
        leaf_names = [leaf_name.split(' ')[0] if not is_ref else leaf_name
                      for leaf_name, is_ref in leaves.items()]
        files = _to_files(leaf_names)
        txt = '\n'.join(files) if as_files else ', '.join(leaf_names)
        mrca = '%d match%s' % (len(leaves), 'es' if len(leaves) > 1 else '')
    else:
        raise ValueError('illegal call')

except (ToytreeError, ValueError) as ex:
    txt = str(ex)
    mrca = ''
    leaves = dict()

# save result in file
_path = path.join('queries', call + '_' + motifs + '_' + ex_motifs)
_path += '_files.txt' if as_files else '_samples.txt'
with open(_path, 'w') as txt_fh:
    txt_fh.write(txt.replace(', ', '\n') + '\n')

# diversity stats
if len(leaves) > 1:

    # check if a reference is among the selected leaves:
    picked_refs = set()
    picked_ids = set()
    for leaf, is_ref in leaves.items():
        if is_ref:
            picked_refs.add(leaf)
        else:
            picked_ids.add(leaf)

    # if a reference is selected
    if picked_refs:
        # read in tabular data
        df = pd.read_csv(form['metadata'].value, sep='\t', dtype={'id': str})
        df.set_index('id', inplace=True)
        df = df.reference_species[pd.notna(df.reference_species)]  # pandas.core.series.Series

        ref_ids = set()
        for k, v in df.items():
            turn_around = ' '.join(v.split(' strain ')[::-1])
            if turn_around in picked_refs:
                picked_refs.remove(turn_around)
                ref_ids.add(k)

        # merge picked_ids and ref_ids
        picked_ids |= ref_ids

    try:
        # per-gene metrics and Tajima's D
        pop_seqs = {record.id.split(' ')[0]: record.seq
                    for record in SeqIO.parse(form['msa_path'].value, 'fasta')
                    if record.id.split(' ')[0] in picked_ids}
        if not pop_seqs:
            raise ValueError('no matches')

        # get gene lengths from form
        g_lens = [entry.split(':') for entry in form['g_lens'].value.split('_')]
        pos = 0
        rows = ''

        # calculate diversity/neutrality statistics per gene
        for entry in g_lens:
            gene, glen = entry[0], int(entry[1])
            cut_seqs = {seq_id: seq[pos:pos + glen] for seq_id, seq in pop_seqs.items()}
            if not cut_seqs:
                raise ValueError('no matches')
            rows += _diversity_stats(gene, cut_seqs, thresholds, poly_allowed)
            pos += glen

        # overall if sensible
        if len(g_lens) > 1:
            rows += _diversity_stats('overall', pop_seqs, thresholds, poly_allowed)

        # write raw html table
        table = '<hr/><h4>Neutrality + Diversity Statistics</h4><p>Exclusion threshold for gap sites ' \
                'was <code>%.2f</code>, and <code>%.2f</code> for unknown sites.</p>' \
                '<table border="1" class="dataframe"><thead><tr style="text-align: justify;"><th></th>' \
                '<th>valid sites</th><th>S</th><th>k</th><th>π</th><th>Watterson\'s θ</th>' \
                '<th>Tajima\'s D</th><th>Genotypes</th><th>gaps</th><th>unknown</th></tr></thead>' \
                '<tbody>%s</tbody></table>' % (thresholds[0], thresholds[1], rows)
    except ValueError:
        table = ''
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
ET.xpath('//link[@id="favicon"]')[0].set('href', 'http://localhost:%d/.img/favicon.ico' % port)
msa_path = ET.xpath('//img[@id="msa_viz"]')
if len(msa_path) > 0:
    msa_path[0].set('src', 'http://localhost:%d/rectangular_msa.png' % port)

# re-set metadata.tsv link
ET.xpath('//a[@id="tsv_link"]')[0].set('href', '../../metadata.tsv')
ET.xpath('//iframe[@id="mview_link"]')[0].set('src', '../../msa_mview.html')
# paste into input fields
ET.xpath('//input[@id="motifs"]')[0].set('value', motifs)
ET.xpath('//input[@id="exclude"]')[0].set('value', ex_motifs)

# create a markup for selected nodes
markup = ';font-weight: bold; font-style: italic; text-decoration: underline 0.15rem'

for g in ET.findall('//svg[@class="toyplot-canvas-Canvas"]/g/g/g[@class="'
                    'toyplot-mark-Text"]/g/g[@class="toyplot-Datum"]/text'):
    st = g.get('style').replace(markup, '')
    # remove old markup
    g.set('style', st)
    if g.text in leaves:
        # insert new markup
        g.set('style', st + markup)

# convert to raw text
print(html.tostring(ET, encoding=sys.stdout.encoding).decode())
