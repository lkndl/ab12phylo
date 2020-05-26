#!/usr/bin/python3
# 2020 Leo Kaindl

"""
Main module of the package. Imports the other package
modules and serves as the interface between them.
"""

import sys
from ab12phylo import cli, its_io, msa, blast, raxml, phylo


def _main():
    """Entry point for primary use"""

    # get command line arguments
    args = cli.parser(sys.argv[1:]).args

    # read in files
    reader = its_io.reader(args)

    # write .FASTAs and source information
    writer = its_io.writer(args, reader)

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

    print('BYE!')


if __name__ == '__main__':
    _main()
