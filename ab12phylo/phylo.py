#!/usr/bin/python3
# 2020 Leo Kaindl

"""
Visualizes tree computed by RAxML-NG using toytree and toyplot.
Renders html using jinja2 and CGI.
"""

import copy
import logging
import os
import pickle
import shutil
import socket
import stat
import subprocess
import sys
import webbrowser
from os import path
from time import time

import numpy as np
import pandas as pd
import random
import toyplot
import toytree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from jinja2 import Template
from lxml import etree
from toyplot import html, png, color, data, locator

from ab12phylo import main, cli

# adapted kxlin colors
raw_rgbas = [(.56, 1, .14, 1), (.16, .44, .8, 1), (.06, 1, .6, 1), (.92, 1, .7, 1),
             (.94, .94, .94, .6), (1, 1, 1, 0), (1, 1, 1, 0), (1, 1, 1, 0)] + [(1, 0, 0, 1)] * 20
kxlin = color.Palette(colors=[color.rgba(i[0], i[1], i[2], i[3]) for i in raw_rgbas])
# author colors
raw_rgbas = [(1.0000, 0.9453, 0.9101), (1.0000, 0.9257, 0.1406), (0.9687, 0.7968, 0.6523),
             (0.5117, 0.4648, 0.6132), (0.1640, 0.6757, 0.9882), (0.3750, 0.4726, 0.7851),
             (0.1171, 0.1718, 0.3242), (0.4960, 0.1484, 0.3281), (1.0000, 0.0039, 0.3007),
             (1.0000, 0.4687, 0.6640)]
pal2 = toyplot.color.Palette(colors=[toyplot.color.rgba(i[0], i[1], i[2], 1) for i in raw_rgbas])

blue = (0.1960, 0.5333, 0.7411, 1)
red = (0.6196, 0.0039, 0.2588, 1)

# seaborn rocket palette
rocket = [blue] + [(0.1237, 0.0717, 0.1822), (0.2452, 0.1049, 0.2639),
                   (0.3809, 0.1206, 0.3250), (0.5172, 0.1179, 0.3545),
                   (0.6582, 0.0955, 0.3536), (0.7965, 0.1050, 0.3106),
                   (0.8949, 0.2178, 0.2531), (0.9429, 0.3754, 0.2636),
                   (0.9592, 0.5330, 0.3748), (0.9644, 0.6702, 0.5150),
                   (0.9689, 0.7980, 0.6851)]
rocket = color.Palette(colors=[color.rgba(t[0], t[1], t[2], 1) for t in rocket])

# seaborn cubehelix palette
cubehelix = [blue] + [[0.1725, 0.1195, 0.2432], [0.2109, 0.1988, 0.3548], [0.2323, 0.2908, 0.4444],
                      [0.2499, 0.3901, 0.5053], [0.2775, 0.4896, 0.5382], [0.3256, 0.5824, 0.5512],
                      [0.3949, 0.6591, 0.5567], [0.4926, 0.7267, 0.5693], [0.6081, 0.7816, 0.6017],
                      [0.7294, 0.8282, 0.6624], [0.8423, 0.8737, 0.7524]]
cubehelix = color.Palette(colors=[color.rgba(t[0], t[1], t[2], 1) for t in cubehelix])

blue = color.rgba(blue[0], blue[1], blue[2], blue[3])
red = color.rgba(red[0], red[1], red[2], red[3])

pg = color.brewer.palette('PinkGreen', reverse=True) \
     + color.Palette(colors=[color.rgba(.26, .26, .26, 1), color.brewer.palette('Spectral')[1]])


def _visualize(*args):
    """Entry point for re-visualization."""
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
        log.info('starting CGI server on port: %d for 10min' % port)
        subprocess.run('python3 -m http.server --cgi %d' % port, shell=True,
                       stdout=subprocess.PIPE, cwd=_dir, timeout=600)
    except subprocess.TimeoutExpired:
        txt = 'CGI server shut down after timeout. If necessary, re-start in %s via ' \
              '"python3 -m http.server --cgi %d" or re-run ab12phylo-view' % (_dir, port)
        log.info(txt)
        print(txt, file=sys.stderr)
    except KeyboardInterrupt:
        txt = 'AB12PHYLO shut down.'
        log.info(txt)
        print('\n' + txt, file=sys.stderr)
    return


class tree_build:
    """
    Tree visualization based on toytree and toyplot. Also calls tree_view()
    """

    def __init__(self, _args):
        self.log = logging.getLogger(__name__)
        self.args = _args
        self.replace = _args.replace

        self.scale = 0.1 if _args.metric == 'FBP' else 10
        self.args.threshold *= 10
        self.tree_file = _args.final_tree + '_%s.nwk' % _args.metric

        preview = self._preview_topology()

        render_dict = self._annotate_msa()

        colors = self._edit1()

        self.log.debug('drawing circular tree')
        start = time()
        ccanvas, axc = self.tree.draw(width=720, height=800, scalebar=True,
                                      node_sizes=colors[2], node_colors=colors[3],
                                      tip_labels=True, tip_labels_align=False, tip_labels_colors=colors[1],
                                      layout='c', edge_type='c')
        ccanvas.style['background-color'] = 'white'
        axc.show = False
        png.render(ccanvas, path.join(self.args.dir, 'circular.png'), scale=1.6)
        self.log.debug('rendered circular in %.2f sec' % (time() - start))

        # prep rectangular tree: save species and pid to name
        tree_no_msa = self._edit2(copy.deepcopy(self.tree))
        # write new-ish newick file
        tree_no_msa.write(self.args.annotated_tree, tree_format=0)

        self.log.debug('drawing tree w/o msa')
        start = time()
        # get dim for canvas
        w, h = 1200, len(self.tips) * 14 + 80
        rcanvas, axes = tree_no_msa.draw(width=w, height=h, scalebar=True, tip_labels=True,
                                         node_sizes=colors[2], node_colors=colors[3],
                                         tip_labels_colors=colors[1])
        rcanvas.style['background-color'] = 'white'
        axes.y.show = False
        axes.x.show = True
        axes.x.domain.max = self.tree.treenode.height / 5  # 0 is right-most tip of tree. -> divide space!

        png.render(rcanvas, path.join(self.args.dir, 'rectangular.png'), scale=1.6)
        self.log.debug('rendered rectangular in %.2f sec' % (time() - start))

        self._mview_msa()

        # render rectangular tree with MSA
        if self.args.msa_viz:
            self.log.debug('drawing tree w/ msa')
            start = time()

            rcanvas_msa = self._with_matrix(colors, kxlin)
            png.render(rcanvas_msa, path.join(self.args.dir, 'rectangular_msa.png'), scale=2)
            self.log.debug('rendered w/ msa in %.2f sec' % (time() - start))

        self._write_html(render_dict, preview, ccanvas, rcanvas, tree_no_msa)
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

        # set dimensions of the canvas
        preview = toyplot.Canvas(width=720, height=360, style={'background-color': 'white'})
        # dissect canvas into two cartesian areas
        ax0 = preview.cartesian(bounds=('5%', '48%', '5%', '95%'))
        ax1 = preview.cartesian(bounds=('52%', '95%', '5%', '95%'))
        style = {'tip_labels': False, 'edge_type': 'c'}

        # call draw with the 'axes' argument to pass it to a specific cartesian area
        self.log.debug('drawing preview')
        self.tree.draw(axes=ax0, layout='r', tip_labels=False)
        self.tree.draw(axes=ax1, layout='c', edge_type='c', tip_labels=False)

        # hide the axes coordinates
        ax0.show = False
        ax1.show = False

        # render to imag
        png.render(preview, self.args.topo)
        self.log.info('rendered preview %s' % self.args.topo)
        return preview

    def _annotate_msa(self):
        """
        Reads genes and no_BLAST from log. Reads in metadata for first gene.
        Writes re-named, re-ordered msa and saves records if necessary.

        :return some information for HTML rendering as dict
        """
        log = open(self.args.log, 'r').read()

        # jinja
        render_info = dict()
        render_info['run_start'] = log[:19]
        render_info['run_args'] = log[50:log.find('\n')]

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

        if '--SKIPPING BLAST--' in log:
            self.args.no_BLAST = True

        # fetch tips labels
        self.tips = list(reversed(self.tree.get_tip_labels()))
        # jinja
        render_info['num_seq'] = str(len(self.tips))
        render_info['genes'] = ', '.join(self.genes)

        # read in tabular data
        self.df = pd.read_csv(self.args.tsv, sep='\t', index_col=1)
        self.df.columns = ['gene'] + list(self.df.columns[1:])
        # drop and order rows
        self.df = self.df[self.df.gene == gene]
        self.df = self.df.loc[self.tips]

        # read in MSA
        records = {record.id: record.seq for record in SeqIO.parse(self.args.msa, 'fasta')}

        # write re-named, re-ordered MSA
        with open(self.args.new_msa, 'w') as msa:
            for tip in self.tips:
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
                    tip = entry.accession
                    des = entry.reference_species

                # get replaced IDs
                if 'replaces' in entry and not pd.isnull(entry.replaces):
                    des += '; replaces ' + entry.replaces

                SeqIO.write(SeqRecord(Seq(seq), id=tip, description=des), msa, 'fasta')
            self.log.debug('wrote updated MSA: %s' % self.args.new_msa)

        # save seqs if necessary
        if self.args.msa_viz:
            self.records = records
        return render_info

    def _with_matrix(self, colors, pal):
        """
        Create a toytree with a matrix representation of the MSA.
        """
        start = time()

        # get lengths of genes
        glens = [(gene, len(seq)) for gene, seq in
                 zip(self.genes, self.records[next(iter(self.records))].split(self.args.sep))]
        self.log.debug('gene lengths: %s' % glens)

        # get ticks for gene ends
        ticks = [item[1] for item in glens]
        for i in range(1, len(ticks)):
            ticks[i] += ticks[i - 1] + 12
        ticks = [p - 12 for p in ticks]  # shift sideways
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
        if self.args.no_BLAST:
            dt = lines.fillna('')
        else:
            lines = lines.fillna('no BLAST hit')
            dt = pd.concat([self.df.loc[:, ['pid']], lines.iloc[:, 0]], axis=1)
        # rename last column
        dt.columns = list(dt.columns)[:-1] + ['species']
        dt = toyplot.data.Table(dt)

        hh, ww = coded_seqs.shape

        # create a canvas
        w = int(ww * 2.2 + hh * 2 + 80)
        h = hh * 16 + 80

        self.log.debug('dim w/ msa %spx by %spx' % (w, h))
        rcanvas = toyplot.Canvas(width=w, height=h, style={'background-color': 'white'})
        axes = rcanvas.cartesian(bounds=(40, 0.26 * w, 40, h - 40))  # xstart xend ystart yend

        self.log.info('drawing tree')
        self.tree.draw(axes=axes, scalebar=True, node_sizes=colors[2], node_colors=colors[3],
                       tip_labels=True, tip_labels_align=True, tip_labels_colors=colors[0],
                       edge_align_style={'stroke-width': .7, 'stroke': 'silver'})
        axes.y.show = False
        axes.x.show = True

        self.log.info('drawing matrix')
        msa = rcanvas.matrix((coded_seqs, color.CategoricalMap(palette=pal)),
                             lshow=False, bshow=True, step=100, bounds=(0.25 * w - 40, w, 0, h),
                             tlocator=locator.Explicit(ticks, [i[0] + '|' for i in glens], format='{:>12}'.format))
        msa.body.gaps.rows[...] = 7

        # add BLAST score indicators
        self.log.debug('adding annotations')
        if not self.args.no_BLAST:
            axa = msa.right.column[-2].cartesian()
            axa.scatterplot(np.zeros(len(self.tips)) + 8, np.arange(len(self.tips)),
                            color=colors[4], size=8, marker='s')
            msa.right.column[-2].width = 30

        # write species column
        msa.right.column[-1].data = dt['species']
        msa.right.column[-1].width = 240
        msa.right.column[-1].align = 'left'
        msa.right.column[-1].lstyle = {'font-weight': 'regular'}

        # colorize species annotation for ref and missing
        shift = msa.shape[0] - len(self.tips) - 2  # line index shift because of 'top'
        for j in range(len(self.tips)):
            msa.cells.cell[j + shift, -1].lstyle = {'fill': list(reversed(colors[0]))[j]}
        self.log.info('rectangular tree setup in %.2f sec' % (time() - start))

        return rcanvas

    def _mview_msa(self):
        """
        Using MView, create an .html visualization for the annotated MSA.
        """
        # hacky re-configuration
        mv_path = path.join(path.abspath(path.dirname(__file__)), 'tools', 'mview-1.67')
        perl_binary = shutil.which('perl')
        lines = open(path.join(mv_path, 'bin', 'mview'), 'r').read().split('\n')
        lines[0] = '#!' + perl_binary
        line = lines[25].split('"')
        line[1] = mv_path
        lines[25] = '"'.join(line)

        # over-write file
        open(path.join(mv_path, 'bin', 'mview'), 'w').write('\n'.join(lines))
        self.log.debug('MView config')

        arg = '%s %s -in fasta -moltype dna -width 80 -conservation on -coloring consensus -threshold 80 ' \
              '-label2 -label4 -label5 -consensus on -con_threshold 80 -html head -css on  -colormap kxlin %s > %s' % \
              (perl_binary, path.join(mv_path, 'bin', 'mview'), self.args.new_msa, self.args.mview_msa)

        p = subprocess.run(arg, shell=True, stdout=subprocess.PIPE)
        self.log.debug('MView run ' + p.stdout.decode('utf-8').strip())
        return

    def _edit1(self):
        """
        type -> reference, normal seq or no BLAST hit
        BLAST -> colors change from green to red in 2% steps
        support -> larger dot means better support.
        :return:
        """

        type_colors, BLAST_colors = dict(), dict()

        # rename reference nodes and add features
        for node in self.tree.treenode.traverse():
            if node.is_leaf():
                entry = self.df.loc[node.name]

                if node.name.startswith('REF_'):
                    node.name = entry.accession
                    type_colors[node.idx] = blue
                    node.add_feature('species', entry.reference_species)
                    BLAST_colors[node.idx] = 0  # this is the position of blue in the palette

                elif 'BLAST_species' in entry and not pd.isnull(entry.BLAST_species):
                    type_colors[node.idx] = color.rgba(.2, .2, .2, 1)
                    node.add_feature('species', entry.BLAST_species)
                    # node.add_feature('score', '%.2f' % entry.pid if entry.pid < 100 else '100')
                    node.add_feature('pid', entry.pid)
                    BLAST_colors[node.idx] = 1 + np.ceil(min(100 - entry.pid, 20) / 2)  # map to the palette

                else:
                    type_colors[node.idx] = rocket[-2] if self.args.no_BLAST else red
                    BLAST_colors[node.idx] = 1 if self.args.no_BLAST else 11
                    # dark gray if no BLASTing, else bad light color

                if self.replace and 'replaces' in entry and not pd.isnull(entry.replaces):
                    node.name += ', ' + entry.replaces

        # convert to lists ordered by increasing index -> matching tips
        label_colors = [type_colors[i] for i in sorted(type_colors)]
        score_colors = [rocket[BLAST_colors[i]] for i in sorted(BLAST_colors)]
        score_colors_rg = [cubehelix[BLAST_colors[i]] for i in sorted(BLAST_colors)]

        # use support values for node size
        min_size = 2
        node_sizes = [min_size + float(siz) * self.scale if siz else 0
                      for siz in self.tree.get_node_values('support', 1, 0)]
        node_colors = [kxlin[1] if siz > self.args.threshold + min_size else pal2[8] for siz in node_sizes]
        return label_colors, score_colors, node_sizes, node_colors, score_colors_rg

    def _edit2(self, _tree):
        """
        Appends species information to node.name for visualization without msa.

        :return: an edited copy of the :param: _tree
        """
        # rename all leaves
        for node in _tree.treenode.traverse():
            if node.is_leaf():
                try:
                    node.name += ' ' + node.species
                    node.name += ' %.2f' % node.pid
                except AttributeError:
                    continue
        return _tree

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
        self.log.debug('rendered embeddable HTMLs in %.2f sec' % (time() - start))

        start = time()
        materials['missing'] = open(self.args.missing_samples, 'r').read()
        materials['msa'] = open(self.args.new_msa, 'r').read()
        materials['newick'] = open(self.args.annotated_tree, 'r').read()
        materials['mview'] = open(self.args.mview_msa, 'r').read()
        self.log.debug('read files in %.2f sec' % (time() - start))

        materials['metric'] = self.args.metric
        materials['msa_path'] = path.abspath(self.args.new_msa)
        materials['threshold'] = int(self.args.threshold * 10)
        materials['min_dist'] = self.args.min_dist
        materials['metadata'] = path.relpath(self.args.tsv, self.args.dir)
        materials['version'] = main.__version__
        materials['author'] = main.__author__

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
        df = pd.read_csv(self.args.tsv, sep='\t', index_col=1)
        df.columns = ['gene'] + list(df.columns[1:])
        # convert box to integer
        df.box = df.box.fillna(0).astype(int)
        # replace missing box with -
        df.box = df.box.replace(0, '-')
        # display only two digits
        pd.options.display.float_format = '{:,.2f}'.format
        if not self.args.no_BLAST:
            # print 100 as exact number
            df.pid = df.pid.replace(100.00, '100')
        # print only file names, not full path
        df.file = df.file.apply(lambda _path: _path.split('/')[-1])
        materials['table'] = df.to_html(na_rep='', justify='justify')

        df2 = pd.read_csv(self.args.bad_seqs, sep='\t')
        # re-order columns
        df2 = df2.reindex(columns=['id', 'problem', 'gene', 'box', 'file'])
        materials['bad_seqs'] = df2.to_html(na_rep='', index=False, justify='justify')
        self.log.debug('rendered tables in %.2f sec' % (time() - start))

        # dump tree and function
        with open(path.join(self.args.dir, 'tree.pickle'), 'wb') as out_file:
            pickle.dump(tree, out_file)
        materials['pickle'] = path.abspath(path.join(self.args.dir, 'tree.pickle'))

        # copy cgi-bin directory to output directory
        shutil.copytree(src=path.join(path.dirname(__file__), 'cgi-bin'),
                        dst=path.join(self.args.dir, 'cgi-bin'), dirs_exist_ok=True)
        # copy logo and icon
        _img = path.join(self.args.dir, '.img')
        os.makedirs(_img, exist_ok=True)
        shutil.copy(src=path.join(self.args.dir, 'cgi-bin', 'favicon.png'),
                    dst=path.join(_img, 'favicon.png'))
        shutil.copy(src=path.join(self.args.dir, 'cgi-bin', 'favicon.ico'),
                    dst=path.join(_img, 'favicon.ico'))

        # make CGI script executable
        cgi_file = path.join(self.args.dir, 'cgi-bin', 'ab12phylo.py')
        os.chmod(cgi_file, os.stat(cgi_file).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        self.log.debug('copied files to output directory')

        start = time()
        self.log.debug('rendering output HTML')
        materials['log'] = open(self.args.log, 'r').read()
        template = Template(open(path.join(path.dirname(__file__), 'template.jinja'), 'r').read())
        html = template.render(materials)
        open(path.join(self.args.dir, 'result.html'), 'w').write(html)
        return
