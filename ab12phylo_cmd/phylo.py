#!/usr/bin/env python3
# 2021 Leo Kaindl

"""
Visualizes tree computed by RAxML-NG using toytree and toyplot.
Renders html using jinja2 and CGI. Visualizes MSA using MView.
Contains both non-primary entry points; -viz and -view.
Allows tree modifications.
"""

import copy
import logging
import math
import multiprocessing
import os
import pickle
import random
import shutil
import socket
import stat
import subprocess
import sys
import webbrowser
import xml.etree.ElementTree as ET
from os import path
from pathlib import Path
from time import time

import numpy as np
import pandas as pd
import toyplot
import toytree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from jinja2 import Template
from lxml import etree
from toyplot import html, pdf, png, svg, color, data, locator

from ab12phylo_cmd import cli
from ab12phylo_cmd.__init__ import __version__, __author__

# adapted kxlin colors
kxlin = [(.56, 1, .14, 1), (.16, .44, .8, 1), (.06, 1, .75, 1), (.92, 1, .7, 1),
         (.94, .94, .94, .6), (1, 1, 1, 0), (1, 1, 1, 0), (1, 1, 1, 0)] + [(1, 0, 0, 1)] * 20
kxlin_pal = color.Palette(colors=[color.rgba(i[0], i[1], i[2], i[3]) for i in kxlin])
# author colors
pal2 = [(1.0000, 0.9453, 0.9101), (1.0000, 0.9257, 0.1406), (0.9687, 0.7968, 0.6523),
        (0.5117, 0.4648, 0.6132), (0.1640, 0.6757, 0.9882), (0.3750, 0.4726, 0.7851),
        (0.1171, 0.1718, 0.3242), (0.4960, 0.1484, 0.3281), (1.0000, 0.0039, 0.3007),
        (1.0000, 0.4687, 0.6640)]
# pal2 = toyplot.color.Palette(colors=[toyplot.color.rgba(i[0], i[1], i[2], 1) for i in pal2])

blue = (0.1960, 0.5333, 0.7411)
red = (0.6196, 0.0039, 0.2588)

# seaborn rocket palette
rocket = [blue] + [(0.1237, 0.0717, 0.1822), (0.2452, 0.1049, 0.2639),
                   (0.3809, 0.1206, 0.3250), (0.5172, 0.1179, 0.3545),
                   (0.6582, 0.0955, 0.3536), (0.7965, 0.1050, 0.3106),
                   (0.8949, 0.2178, 0.2531), (0.9429, 0.3754, 0.2636),
                   (0.9592, 0.5330, 0.3748), (0.9644, 0.6702, 0.5150),
                   (0.9689, 0.7980, 0.6851)]
# rocket = color.Palette(colors=[color.rgba(t[0], t[1], t[2], 1) for t in rocket])

# seaborn cubehelix palette
cubehelix = [blue] + [[0.1725, 0.1195, 0.2432], [0.2109, 0.1988, 0.3548], [0.2323, 0.2908, 0.4444],
                      [0.2499, 0.3901, 0.5053], [0.2775, 0.4896, 0.5382], [0.3256, 0.5824, 0.5512],
                      [0.3949, 0.6591, 0.5567], [0.4926, 0.7267, 0.5693], [0.6081, 0.7816, 0.6017],
                      [0.7294, 0.8282, 0.6624], [0.8423, 0.8737, 0.7524]]
# cubehelix = color.Palette(colors=[color.rgba(t[0], t[1], t[2], 1) for t in cubehelix])

# blue = color.rgb(blue[0], blue[1], blue[2])
# red = color.rgb(red[0], red[1], red[2])

pg = color.brewer.palette('PinkGreen', reverse=True) \
     + color.Palette(colors=[color.rgb(.26, .26, .26), color.brewer.palette('Spectral')[1]])

base_legend = '<svg class="toyplot-canvas-Canvas" xmlns:toyplot="http://www.sandia.gov/toyplot" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" width="100.0px" height="40.0px" viewBox="0 0 100.0 40.0" preserveAspectRatio="xMidYMid meet" style="background-color:transparent;border-color:#292724;border-style:none;border-width:1.0;fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:Helvetica;font-size:12px;opacity:1.0;stroke:rgb(16.1%,15.3%,14.1%);stroke-opacity:1.0;stroke-width:1.0" id="ta084f8acf92540e38738f50593e76861"><g class="toyplot-coordinates-Table" id="t9e2f7ff350ea4dbebe0c12a7395bb879"><rect x="0.0" y="0.0" width="20.0" height="20.0" style="fill:rgb(56%,100%,14%);fill-opacity:1.0;stroke:rgb(100%,100%,100%);stroke-opacity:1.0;stroke-width:4" /><rect x="20.0" y="0.0" width="20.0" height="20.0" style="fill:rgb(16%,44%,80%);fill-opacity:1.0;stroke:rgb(100%,100%,100%);stroke-opacity:1.0;stroke-width:4" /><rect x="40.0" y="0.0" width="20.0" height="20.0" style="fill:rgb(6%,100%,75%);fill-opacity:1.0;stroke:rgb(100%,100%,100%);stroke-opacity:1.0;stroke-width:4" /><rect x="60.0" y="0.0" width="20.0" height="20.0" style="fill:rgb(92%,100%,70%);fill-opacity:1.0;stroke:rgb(100%,100%,100%);stroke-opacity:1.0;stroke-width:4" /><rect x="80.0" y="0.0" width="20.0" height="20.0" style="fill:rgb(94%,94%,94%);fill-opacity:0.6;stroke:rgb(100%,100%,100%);stroke-opacity:1.0;stroke-width:4" /><g transform="translate(10.0,30.0)"><text x="-4.332" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">A</text></g><g transform="translate(30.0,30.0)"><text x="-4.332" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">C</text></g><g transform="translate(50.0,30.0)"><text x="-4.668" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">G</text></g><g transform="translate(70.0,30.0)"><text x="-3.666" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">T</text></g><g transform="translate(90.0,30.0)"><text x="-4.332" y="3.066" style="fill:rgb(16.1%,15.3%,14.1%);fill-opacity:1.0;font-family:helvetica;font-size:12.0px;font-weight:bold;stroke:none;vertical-align:baseline;white-space:pre">N</text></g></g></svg>'


def _visualize(*args):
    """Entry point for re-visualization"""
    # if no argument -> default to cwd
    args = sys.argv[1:]
    if len(args) == 0:
        args = [os.getcwd()]
    args = ['--visualize'] + args
    # parse to get defaults
    namespace = cli.parser(args, visualize=True).args
    tree_build(namespace)
    print('BYE!')


def _view():
    """Entry point for re-view"""
    # if no argument -> default to cwd
    args = sys.argv[1:]
    if len(args) == 0:
        args = [os.getcwd()]
    args = ['-view'] + args
    # parse to get defaults
    namespace = cli.parser(args, view=True).args
    tree_view(namespace.dir)
    print('BYE!')


def tree_view(_dir):
    """
    Opens a result page and starts a CGI server
    :param _dir:
    :return:
    """
    log = logging.getLogger(__name__)

    # open in new tab
    webbrowser.open(path.join(_dir, 'result.html'), new=2)

    # find original port
    port = int(etree.parse(path.join(_dir, 'result.html'), parser=etree.HTMLParser())
               .findall('//meta[@name="port"]')[0].get('content'))

    # start CGI server
    try:
        py3 = sys.executable
        py3 = py3 if ' ' not in py3 else f'"{py3}"'
        log.debug('python3 executable at %s' % py3)
        log.info('starting CGI server on port: %d for 20min' % port)
        subprocess.run('%s -m http.server --cgi %d' % (py3, port), shell=True,
                       stdout=subprocess.PIPE, cwd=_dir, timeout=1200)
    except subprocess.TimeoutExpired:
        txt = 'CGI server shut down after timeout. If necessary, re-start in %s via ' \
              '"%s -m http.server --cgi %d" or re-run ab12phylo-view' % (_dir, sys.executable, port)
        log.info(txt)
        print(txt, file=sys.stderr)
    except KeyboardInterrupt:
        txt = 'AB12PHYLO shut down.'
        log.info(txt)
        print('\n' + txt, file=sys.stderr)
    return


def mview_msa(args):
    """
    Using MView, create an .html visualization for the annotated MSA.
    """
    # hacky re-configuration
    mv_path = path.join(path.abspath(path.dirname(__file__)),
                        'tools', 'mview-1.67').replace(' ', '\ ')
    perl_binary = shutil.which('perl')
    if not perl_binary:
        raise ValueError
    lines = open(path.join(mv_path, 'bin', 'mview'), 'r').read().split('\n')
    # monkey patch
    lines[0] = '#!' + perl_binary
    line = lines[25].split('"')
    line[1] = mv_path
    lines[25] = '"'.join(line)

    # over-write file
    open(path.join(mv_path, 'bin', 'mview'), 'w').write('\n'.join(lines))

    used_msa = args.new_msa if path.isfile(args.new_msa) else args.msa

    arg = '%s %s -in fasta -moltype dna -width 80 -conservation on -coloring consensus ' \
          '-threshold 80 -label2 -label4 -label5 -consensus on -con_threshold 80 ' \
          '-html head -css on  -colormap leo "%s" > "%s"' % \
          (perl_binary, path.join(mv_path, 'bin', 'mview'), used_msa, args.mview_msa)

    p = subprocess.run(arg, shell=True, stdout=subprocess.PIPE)
    return 'MView run ' + p.stdout.decode('utf-8').strip()


class tree_build:
    """
    Tree visualization based on toytree and toyplot. Also calls tree_view()
    """

    def __init__(self, _args):
        self.log = logging.getLogger(__name__)
        self.args = _args
        self.replace = _args.replace

        self.args.threshold *= 10
        self.tree_file = _args.final_tree + '_%s.nwk' % _args.metric

        self.log.debug('toyplot %s' % toyplot.__version__)
        self.log.debug('toytree %s' % toytree.__version__)
        self.log.debug(sys.version_info)

        preview = self._preview_topology()

        render_dict = self._annotate_msa()

        self._edit1()

        self.log.debug('drawing circular tree')
        start = time()
        try:
            ctup = self.tree.draw(width=800, height=800, scalebar=True,
                                  node_sizes=list(self.tree.get_node_values('size', 1, 1)),
                                  node_colors=[color.rgb(n[0], n[1], n[2]) for n in list(
                                      self.tree.get_node_values('color', 1, 1))],
                                  tip_labels=True, tip_labels_align=False,
                                  tip_labels_colors=[color.rgb(rocket[n][0], rocket[n][1], rocket[n][2])
                                                     for n in list(self.tree.get_node_values('score', 1, 1))
                                                     if n != -1][::-1],
                                  layout='c')
            ctup[0].style['background-color'] = 'white'
            ctup[1].show = False
            if not self.args.out_fmt:
                self.args.out_fmt = ['pdf']

            if 'png' in self.args.out_fmt:
                png.render(ctup[0], path.join(self.args.dir, 'circular.png'), scale=1.6)
            if 'pdf' in self.args.out_fmt or len(self.args.out_fmt) == 0:
                pdf.render(ctup[0], path.join(self.args.dir, 'circular.pdf'))
            if 'svg' in self.args.out_fmt:
                svg.render(ctup[0], path.join(self.args.dir, 'circular.svg'))
            self.log.info('rendered circular in %.2f sec' % (time() - start))

        except ET.ParseError as ex:
            self.log.error('XML ParseError. Invalid character in a sample ID? Please check metadata.tsv')
            sys.exit(1)
        except ValueError as ex:
            self.log.exception(ex)

        # prep rectangular tree: save species and pid to name
        tree_no_msa = copy.deepcopy(self.tree)

        # save species and pid to name
        for node in tree_no_msa.treenode.traverse():
            if node.is_leaf():
                try:
                    if node.name in node.species:
                        node.name = node.species
                    else:
                        node.name += ' ' + node.species
                        node.name += ' %.2f' % node.pid
                except AttributeError:
                    continue

        # write new-ish newick file
        tree_no_msa.write(self.args.annotated_tree, tree_format=0)

        self.log.debug('drawing tree without msa')
        start = time()
        try:
            # get dim for canvas
            w, h = 1200, len(self.tips) * 14 + 80
            if self.args.print_supports:
                rtup = tree_no_msa.draw(width=w, height=h, scalebar=True, tip_labels=True, node_labels='support',
                                        node_labels_style={'font-size': '6px', 'fill': '#FFFFFF',
                                                           'baseline-shift': '-1px',
                                                           'font-weight': 'bold'},
                                        node_sizes=list(self.tree.get_node_values('size', 1, 1)),
                                        node_colors=[color.rgb(n[0], n[1], n[2]) for n in
                                                     list(self.tree.get_node_values('color', 1, 1))],
                                        tip_labels_colors=[color.rgb(rocket[n][0], rocket[n][1], rocket[n][2])
                                                           for n in list(self.tree.get_node_values('score', 1, 1))
                                                           if n != -1][::-1])
            else:
                rtup = tree_no_msa.draw(width=w, height=h, scalebar=True, tip_labels=True,
                                        node_sizes=list(self.tree.get_node_values('size', 1, 1)),
                                        node_colors=[color.rgb(n[0], n[1], n[2]) for n in
                                                     list(self.tree.get_node_values('color', 1, 1))],
                                        tip_labels_colors=[color.rgb(rocket[n][0], rocket[n][1], rocket[n][2])
                                                           for n in list(self.tree.get_node_values('score', 1, 1))
                                                           if n != -1][::-1])
            rtup[0].style['background-color'] = 'white'
            rtup[1].y.show = False
            rtup[1].x.show = True
            rtup[1].x.domain.max = self.tree.treenode.height / 5  # 0 is right-most tip of tree. -> divide space!

            if 'png' in self.args.out_fmt:
                png.render(rtup[0], path.join(self.args.dir, 'rectangular.png'), scale=1.6)
            if 'pdf' in self.args.out_fmt:
                pdf.render(rtup[0], path.join(self.args.dir, 'rectangular.pdf'))
            if 'svg' in self.args.out_fmt:
                svg.render(rtup[0], path.join(self.args.dir, 'rectangular.svg'))
            self.log.info('rendered rectangular in %.2f sec' % (time() - start))

        except ValueError as ex:
            self.log.exception(ex)

        try:
            self.log.info(mview_msa(self.args))
        except Exception as ex:
            with open(self.args.mview_msa, 'w') as fh:
                fh.write('MView (Perl) failed\n')
                fh.write(str(ex))

        # render rectangular tree with MSA
        reserve_gap = False
        if type(self.args.msa_viz) == list and 'png' in self.args.msa_viz:
            try:
                reserve_gap = True
                proc = multiprocessing.Process(target=self._with_matrix,
                                               args=([kxlin_pal]))
                proc.start()
            except KeyboardInterrupt:
                self.log.warning('cancel msa_viz')

        # but display it any way if it exists
        if reserve_gap or path.isfile(path.join(self.args.dir, 'rectangular_msa.png')):
            render_dict['msa_viz'] = '<div id="msa_viz" href="msa_viz" class="section level3">' \
                                     '<h3><strong>3.2</strong> +MSA</h3>' \
                                     '<img id="msa_viz" src="rectangular_msa.png" alt="tree+MSA" ' \
                                     'style="width:860px"/><h4>%s</h4></div>' % base_legend
        else:
            render_dict['msa_viz'] = ''

        self._write_html(render_dict, preview, ctup[0], rtup[0], tree_no_msa)
        if not self.args.headless:
            tree_view(self.args.dir)

    def _preview_topology(self):
        """
        Reads in the tree file, and saves a .png with rectangular and circular topology previews.

        :return: toyplot.Canvas object
        """
        self.log.debug('reading tree')
        # read in the tree
        self.tree = toytree.tree(open(self.tree_file, 'r').read(), tree_format=0)

        # drop / replace nodes
        if self.args.replace_nodes:
            for idx in self.args.replace_nodes:
                self.tree.idx_dict[idx].add_sister(name='%d_replaced' % idx, dist=1)
                self.tree.idx_dict[idx].detach()
            # saving the node names now saves some work
            self.args.replace_nodes = ['%d_replaced' % idx for idx in self.args.replace_nodes]
        else:
            # save empty list
            self.args.replace_nodes = []
        if self.args.drop_nodes:
            [self.tree.idx_dict[idx].detach() for idx in self.args.drop_nodes]

        # outgroup rooting
        if self.args.root:
            self.tree = self.tree.root(names=[self.tree.idx_dict[self.args.root].name])

        self.tree = toytree.tree(self.tree.write())

        # set dimensions of the canvas
        preview = toyplot.Canvas(width=800, height=400, style={'background-color': 'white'})
        # dissect canvas into two cartesian areas
        ax0 = preview.cartesian(bounds=('5%', '48%', '5%', '95%'))
        ax1 = preview.cartesian(bounds=('52%', '95%', '5%', '95%'))

        # call draw with the 'axes' argument to pass it to a specific cartesian area
        self.log.debug('drawing preview')
        self.tree.draw(axes=ax0, layout='r', tip_labels=False)
        self.tree.draw(axes=ax1, layout='c', tip_labels=False)

        # hide the axes coordinates
        ax0.show = False
        ax1.show = False

        # render to image
        png.render(preview, self.args.topo)
        self.log.info('rendered preview %s' % self.args.topo)
        return preview

    def _annotate_msa(self):
        """
        Reads genes and no_BLAST from log. Reads in metadata for first gene.
        Writes re-named, re-ordered msa and saves records if necessary.

        :return some information for HTML rendering as dict
        """
        if self.args.finish:
            try:
                log = open(self.args.log[:-5] + '1.log', 'r').read()
                self.args.log = self.args.log[:-5] + '1.log'
            except FileNotFoundError:
                try:
                    log = open(self.args.log[:-7] + '.log', 'r').read()
                    self.args.log = self.args.log[:-7] + '.log'
                except FileNotFoundError:
                    self.log.error('no log file found')
                    exit(1)
        else:
            try:
                log = open(self.args.log, 'r').read()
            except FileNotFoundError:
                log = open(self.args.log[:-4] + '-p1.log', 'r').read()
                self.args.log = self.args.log[:-4] + '-p1.log'

        # jinja
        render_info = dict()
        render_info['run_start'] = log[:19]
        render_info['run_args'] = log[54:log.find('\n')]

        # extract genes
        start = log.find('--GENES--') + 10
        end = log[start:start + 800].find('\n')
        self.genes = log[start:start + end].split('::')
        self.log.debug('read genes from log: %s' % ':'.join(self.genes))
        gene = self.genes[0]
        self.log.debug('using annotation for %s' % gene)

        # extract seed
        start = log.find('seed for this run:') + 19
        end = log[start:start + 20].find('\n')
        render_info['seed'] = int(log[start:start + end])

        if '--SKIPPING BLAST--' in log or 'BLAST+ not installed' in log:
            self.args.no_BLAST = True

        # fetch tips labels
        self.tips = list(reversed(self.tree.get_tip_labels()))
        # jinja
        render_info['num_seq'] = str(len(self.tips))
        render_info['genes'] = ', '.join(self.genes)

        # read in tabular data
        self.df = pd.read_csv(self.args.tsv, sep='\t', dtype={'id': str})
        self.df.set_index('id', inplace=True)
        # drop and order rows
        self.df = self.df[self.df.gene == gene]
        # add dummy data to data frame
        for node_id in self.args.replace_nodes:
            self.df.at[node_id, 'gene'] = gene
        self.df = self.df.reindex(self.tips)

        # read in MSA
        records = {record.id: record.seq for record in SeqIO.parse(self.args.msa, 'fasta')}

        # write re-named, re-ordered MSA
        with open(self.args.new_msa, 'w') as new_msa:
            for tip in self.tips:
                if tip in self.args.replace_nodes:
                    continue
                # get seq and cut out artificial separator
                seq = str(records[tip]).replace(self.args.sep, '')
                # get metadata
                entry = self.df.loc[tip]

                des = ''
                # get BLAST info
                if 'BLAST_species' in entry and not pd.isnull(entry.BLAST_species):
                    des = entry.BLAST_species
                    # annotate BLAST score
                    pid = '%.2f' % entry.pid if entry.pid < 100 else '100'
                    des += '; pid %s' % pid

                # rename REFs
                if tip.startswith('REF_'):
                    if len(self.genes) > 1 and 'strain' in entry.reference_species:
                        des, tip = entry.reference_species.split(' strain ')
                    else:
                        tip = entry.accession
                        des = entry.reference_species

                # get replaced IDs
                if 'replaces' in entry and not pd.isnull(entry.replaces):
                    des += '; replaces ' + entry.replaces

                SeqIO.write(SeqRecord(Seq(seq), id=tip, description=des), new_msa, 'fasta')
            self.log.debug('wrote updated MSA: %s' % self.args.new_msa)

        # get lengths of genes
        self.g_lens = [(gene, len(seq)) for gene, seq in
                       zip(self.genes, records[next(iter(records))].split(self.args.sep))]

        # save seqs if necessary
        if type(self.args.msa_viz) == list:
            self.records = records
            # make empty dummy seqs
            seq = '-' * [item[1] for item in self.g_lens if item[0] == gene][0]
            for _id in self.args.replace_nodes:
                self.records[_id] = Seq(seq)

        self.log.debug('gene lengths: %s' % self.g_lens)
        return render_info

    def _with_matrix(self, pal):
        """
        Create a toytree with a matrix representation of the MSA.
        """
        start = time()
        self.log.warning('drawing tree with msa. You can interrupt and proceed via Cmd+C.')
        start = time()

        # if no format was specified, render .PNGs
        if not self.args.msa_viz:
            self.args.msa_viz = ['png']

        # get ticks for gene ends
        ticks = [item[1] for item in self.g_lens]
        for i in range(1, len(ticks)):
            ticks[i] += ticks[i - 1] + 10
        # ticks[-1] -= 1
        ticks = [p - 3 for p in ticks]  # shift sideways
        self.log.debug('ticks: %s' % ':'.join(str(tick) for tick in ticks))

        # get seqs as numpy array in plot order
        seqs = np.array([self.records[_id] for _id in self.tips])
        self.log.info('numpy array representation of MSA\n\t'
                      'shape: %s entries: %s' % (seqs.shape, ':'.join(np.unique(seqs))))
        bases = ['A', 'C', 'G', 'T', 'N', '-', ' ', 'S']
        weird_chars = np.setdiff1d(np.unique(seqs), bases)
        self.log.debug('observed: %s' % ':'.join(set(list(weird_chars) + bases)))
        # get a mapping from characters to integer
        codes = dict(zip(list(bases) + list(weird_chars), range(len(bases) + len(weird_chars))))
        # translate seqs array to integers
        coded_seqs = np.vectorize(codes.get)(seqs)

        # get a data table of other info
        cols = set(self.df.columns).intersection({'reference_species', 'BLAST_species'})
        lines = self.df[cols]

        if lines.shape[1] == 2:
            # re-order columns
            lines = lines.reindex(columns=['BLAST_species', 'reference_species'])
            # paste ref species for references
            for j in range(lines.shape[0]):
                if not pd.isna(lines.iloc[j, 1]):
                    lines.iloc[j, 0] = lines.iloc[j, 1]

        # annotate missing species
        if self.args.no_BLAST or 'pid' not in self.df:
            dt = lines.fillna('')
        else:
            lines = lines.fillna('no BLAST hit')
            dt = pd.concat([self.df.loc[:, ['pid']], lines.iloc[:, 0]], axis=1)
        # rename last column
        try:
            dt.columns = list(dt.columns)[:-1] + ['species']
            dt = toyplot.data.Table(dt)
        except ValueError:
            dt['species'] = ''

        hh, ww = coded_seqs.shape

        # create a canvas
        w = int(ww * 2.2 + hh * 2 + 80)
        h = hh * 16 + 80

        self.log.debug('dim w/ msa %spx by %spx' % (w, h))
        rcanvas = toyplot.Canvas(width=w, height=h, style={'background-color': 'white'})
        axes = rcanvas.cartesian(bounds=(40, 0.26 * w, 40, h - 40))  # xstart xend ystart yend

        self.log.debug('drawing tree')
        if self.args.print_supports:
            self.tree.draw(axes=axes, scalebar=True, tip_labels=True, tip_labels_align=True,
                           node_labels='support', node_labels_style={'font-size': '6px',
                                                                     'fill': '#FFFFFF',
                                                                     'baseline-shift': '-1px',
                                                                     'font-weight': 'bold'},
                           edge_align_style={'stroke-width': .7, 'stroke': 'silver'},
                           node_sizes=list(self.tree.get_node_values('size', 1, 1)),
                           node_colors=[color.rgb(n[0], n[1], n[2]) for n in list(
                               self.tree.get_node_values('color', 1, 1))],
                           tip_labels_colors=[color.rgb(n[0], n[1], n[2]) for n in list(
                               self.tree.get_node_values('type', 1, 1)) if n != 0][::-1])
        else:
            self.tree.draw(axes=axes, scalebar=True, tip_labels=True, tip_labels_align=True,
                           edge_align_style={'stroke-width': .7, 'stroke': 'silver'},
                           node_sizes=list(self.tree.get_node_values('size', 1, 1)),
                           node_colors=[color.rgb(n[0], n[1], n[2]) for n in list(
                               self.tree.get_node_values('color', 1, 1))],
                           tip_labels_colors=[color.rgb(n[0], n[1], n[2]) for n in list(
                               self.tree.get_node_values('type', 1, 1)) if n != 0][::-1])
        axes.y.show = False
        axes.x.show = True

        self.log.info('drawing matrix')
        msa = rcanvas.matrix((coded_seqs, color.CategoricalMap(palette=pal)),
                             lshow=False, bshow=True, step=100, bounds=(0.26 * w - 40, w, 0, h),
                             tlocator=locator.Explicit(ticks, [i[0] + ']' for i in self.g_lens],
                                                       format='{:>12}'.format))
        msa.body.gaps.rows[...] = 7

        # # add BLAST score indicators
        # self.log.debug('adding annotations')
        # if not self.args.no_BLAST:
        #     axa = msa.right.column[-2].cartesian()
        #     axa.scatterplot(np.zeros(len(self.tips)) + 8, np.arange(len(self.tips)),
        #                     color=colors[4], size=8, marker='s')
        #     msa.right.column[-2].width = 30
        msa.right.column[-2].width = 12
        msa.left.column[1].width = 44

        # write species column
        msa.right.column[-1].data = dt['species']
        msa.right.column[-1].width = 240
        msa.right.column[-1].align = 'left'
        msa.right.column[-1].lstyle = {'font-weight': 'regular'}

        # colorize species annotation for ref and missing
        shift = msa.shape[0] - len(self.tips) - 2  # line index shift because of 'top'
        for j in range(len(self.tips)):
            msa.cells.cell[j + shift, -1].lstyle = \
                {'fill': [color.rgb(cubehelix[n][0], cubehelix[n][1], cubehelix[n][2])
                          for n in list(self.tree.get_node_values('score', 1, 1))
                          if n != -1][j]}
        self.log.info('rectangular tree setup in %.2f sec' % (time() - start))

        if 'png' in self.args.msa_viz:
            png.render(rcanvas, path.join(self.args.dir, 'rectangular_msa.png'), scale=2)
        if 'pdf' in self.args.msa_viz:
            pdf.render(rcanvas, path.join(self.args.dir, 'rectangular_msa.pdf'))
        self.log.info('rendered with msa in %.2f sec' % (time() - start))

    def _edit1(self):
        """
        type -> reference, normal seq or no BLAST hit
        BLAST -> colors change from black to red to glowy yellow in 2% steps
        support -> larger dot means better support.
        :return:
        """
        multi = True if len(self.genes) > 1 else False

        # guess if support val.s are in [0-1] or in [0-100], then scale with 1 or 1/100 resp.
        sup_scaler = 1 if 1.1 > max([float(s) for s in
                                    self.tree.get_node_values('support') if s]) else 100
        # rename reference nodes and add features
        for node in self.tree.treenode.traverse():
            if node.is_leaf():
                entry = self.df.loc[node.name]

                if node.name.startswith('REF_'):
                    if multi and 'strain' in entry.reference_species:
                        entry.at['reference_species'], node.name = entry.reference_species.split(' strain ')
                    else:
                        node.name = entry.accession
                    node.add_feature('type', blue)
                    node.add_feature('species', entry.reference_species)
                    node.add_feature('score', 0)  # this is the position of blue in the palette

                elif 'BLAST_species' in entry and not pd.isnull(entry.BLAST_species):
                    node.add_feature('type', (.2, .2, .2, 1))
                    node.add_feature('species', entry.BLAST_species)
                    # node.add_feature('score', '%.2f' % entry.pid if entry.pid < 100 else '100')
                    node.add_feature('pid', entry.pid)
                    node.add_feature('score', 1 + math.ceil(min(100 - entry.pid, 20) / 2))  # map to the palette

                else:
                    node.add_feature('type', rocket[1] if self.args.no_BLAST else red)
                    node.add_feature('score', 1 if self.args.no_BLAST else 11)
                    # dark gray if no BLASTing, else bad light color

                if self.replace and 'replaces' in entry and not pd.isnull(entry.replaces):
                    node.name += ', ' + entry.replaces

                node.add_feature('size', 0)
                node.add_feature('color', kxlin[1])  # just so the field exists
            else:
                # use support values for node size
                min_size = 2
                node.add_feature('size', min_size + float(node.support) / sup_scaler * 10 if node.support else 0)
                node.add_feature('color', kxlin[1] if node.size > self.args.threshold + min_size else pal2[8])
                node.add_feature('score', -1)  # just so the field exists
                node.add_feature('type', 0)  # just so the field exists

            # some tree tweaks
            if self.args.print_supports:
                node.support = int(round(node.support * 100 / sup_scaler))
            if self.args.min_plot_dist:
                node.dist = max(node.dist, self.args.min_plot_dist)
        return

    def _write_html(self, materials, topo, circ, rect, tree):
        """
        Write interactive output .HTML using Jinja2
        :return:
        """
        self.log.debug('creating HTML')

        start = time()
        materials['topo'] = toyplot.html.tostring(topo)
        materials['circular'] = toyplot.html.tostring(circ)
        materials['rectangular'] = toyplot.html.tostring(rect)
        self.log.info('rendered embeddable HTMLs in %.2f sec' % (time() - start))

        start = time()
        materials['missing'] = open(self.args.missing_samples, 'r').read() \
            if os.path.exists(self.args.missing_samples) else 'no missing samples'
        materials['msa'] = open(self.args.new_msa, 'r').read()
        materials['newick'] = open(self.args.annotated_tree, 'r').read()
        # fully insert MView HTML
        # materials['mview'] = open(self.args.mview_msa, 'r').read()
        self.log.debug('read files in %.2f sec' % (time() - start))

        materials['metric'] = self.args.metric
        materials['gap'] = self.args.gap_share
        materials['unknown'] = self.args.unknown_share
        materials['poly'] = self.args.poly_allelic
        materials['g_lens'] = '_'.join([entry[0] + ':' + str(entry[1]) for entry in self.g_lens])
        materials['msa_path'] = path.basename(self.args.msa)  # path.abspath(self.args.new_msa)
        materials['mview'] = path.basename(self.args.mview_msa)
        materials['threshold'] = int(self.args.threshold * 10)
        materials['min_dist'] = self.args.min_dist
        materials['metadata'] = path.relpath(self.args.tsv, self.args.dir)
        materials['version'] = __version__
        materials['author'] = __author__

        # find a free port. susceptible to race conditions
        sock = socket.socket()
        try:
            sock.bind(('localhost', random.randint(8000, 8010)))
        except OSError:
            sock.bind(('localhost', 0))
        port = sock.getsockname()[1]
        sock.close()
        materials['port'] = port

        start = time()
        df = pd.read_csv(self.args.tsv, sep='\t', dtype={'id': str})
        df.set_index('id', inplace=True)
        # convert box to integer
        df.box = df.box.fillna(0).astype(int)
        # replace missing box with -
        df.box = df.box.replace(0, '-')
        # display only two digits
        pd.options.display.float_format = '{:,.2f}'.format
        if 'pid' in df:
            # print 100 as exact number
            df.pid = df.pid.replace(100.00, '100')
        # print only file names, not full path
        df.file = df.file.fillna('')
        df.file = df.file.apply(lambda _path: _path.split('/')[-1])
        materials['table'] = df.to_html(na_rep='', justify='justify')

        df2 = pd.read_csv(self.args.bad_seqs, sep='\t', dtype={'id': str})
        # re-order columns
        df2 = df2.reindex(columns=['id', 'problem', 'gene', 'box', 'file'])
        df.box = df.box.replace(0, '-')
        materials['bad_seqs'] = df2.to_html(na_rep='', index=False, justify='justify')
        self.log.debug('rendered tables in %.2f sec' % (time() - start))

        # dump tree
        with open(path.join(self.args.dir, 'tree.pickle'), 'wb') as out_file:
            pickle.dump(tree, out_file)
        materials['pickle'] = 'tree.pickle'  # path.abspath(path.join(self.args.dir, 'tree.pickle'))

        # copy cgi-bin directory to output directory
        _dst = path.join(self.args.dir, 'cgi-bin')
        # to keep compatibility with python3.6, don't use shutil.copytree(dirs_exist_ok)
        try:
            # delete old directory
            shutil.rmtree(path=_dst)
        except FileNotFoundError:
            pass
        shutil.copytree(src=path.join(path.dirname(__file__), 'cgi-bin'),
                        dst=path.join(self.args.dir, 'cgi-bin'))

        # copy logo and icon
        _img = path.join(self.args.dir, '.img')
        os.makedirs(_img, exist_ok=True)
        shutil.copy(src=path.join(_dst, 'favicon.png'), dst=path.join(_img, 'favicon.png'))
        shutil.copy(src=path.join(_dst, 'favicon.ico'), dst=path.join(_img, 'favicon.ico'))
        self.log.debug('copied files to output directory')

        start = time()
        self.log.debug('rendering output HTML')
        materials['log'] = open(self.args.log, 'r').read()
        template = Template(open(path.join(path.dirname(__file__), 'template.jinja'), 'r').read())
        rendered_html = template.render(materials)

        # make CGI script executable a bit later so it's surely there
        cgi_file = Path(self.args.dir) / 'cgi-bin' / 'ab12phylo.py'
        cgi_file.chmod(cgi_file.stat().st_mode | stat.S_IEXEC | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        open(path.join(self.args.dir, 'result.html'), 'w').write(rendered_html)
        return
