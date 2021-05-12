#!/usr/bin/env python3
# 2021 Leo Kaindl

"""
Main module of the pipeline. Imports the other package modules and serves as the interface
between them. The `--visualize` and `--view` entry points are in the :class:`phylo` module.
"""

import sys

from ab12phylo_cmd import cli, i_o, msa, blast, ml, phylo


def _main():
    # get command line arguments
    args = cli.parser(sys.argv[1:]).args

    if not args.finish:
        if args.add_xml:
            print('starting -px BLAST XML run', file=sys.stderr)
            blaster = blast.blast_build(args, None)
            blaster.start()
            # reading an XML requires no waiting
            blaster.join()
            # visualize MSA for inspection
            phylo.mview_msa(args)
            print('finished -px BLAST XML run', file=sys.stderr)
            exit(0)

        if args.prepare:
            print('starting -p1 prep run', file=sys.stderr)

        # read in files
        reader = i_o.reader(args)

        # write .FASTAs and source information
        writer = i_o.writer(args, reader)

        # build & trim MSAs
        aligner = msa.msa_build(args, reader.seq_counts)

        # run BLAST and parse results de-synced
        blaster = blast.blast_build(args, reader)
        blaster.start()

        if args.prepare:
            # wait for BLAST if necessary
            blaster.join()
            # visualize MSA for inspection
            phylo.mview_msa(args)
            print('finished -p1 prep run', file=sys.stderr)
            exit(0)

        # build trees with raxml-ng or iqtree2
        ml_threads = ml.ml_build(args)
        ml_threads.run()

        # wait for BLAST if necessary
        blaster.join()

        # visualize best tree
        phylo.tree_build(args)

    else:
        print('starting -p2 finishing run', file=sys.stderr)
        # build trees with raxml-ng or iqtree2
        ml_threads = ml.ml_build(args)
        ml_threads.run()

        # visualize best tree
        phylo.tree_build(args)


if __name__ == '__main__':
    _main()
