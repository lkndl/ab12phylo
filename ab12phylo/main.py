#!/usr/bin/python3
# 2020 Leo Kaindl

"""
Main module of the pipeline. Imports the other package modules and serves as the interface
between them. The `--visualize` and `--view` entry points are in the :class:`phylo` module.
"""
__author__ = 'Leo Kaindl'
__email__ = 'leo.kaindl@tum.de'
__version__ = '0.1b.11'
__date__ = '9 July 2020'
__license__ = 'MIT'
__status__ = 'Beta'

import sys
from ab12phylo import cli, i_o, msa, blast, raxml, phylo


def _main():
    # get command line arguments
    args = cli.parser(sys.argv[1:]).args

    # read in files
    reader = i_o.reader(args)

    # write .FASTAs and source information
    writer = i_o.writer(args, reader)

    # build & trim MSAs
    aligner = msa.msa_build(args, reader.seq_counts)

    # run BLAST and parse results de-synced
    blaster = blast.blast_build(args, reader)
    blaster.start()

    # build trees with raxml-ng
    raxml_threads = raxml.raxml_build(args, reader.metadata)
    raxml_threads.run()

    # wait for BLAST if necessary
    blaster.join()

    # visualize best tree
    phylo.tree_build(args)


if __name__ == '__main__':
    _main()
