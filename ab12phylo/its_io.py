# 2020 Leo Kaindl

"""
This module reads sequencing data to memory using the :class:`its_io.reader`,
does some quality control using :class:`filter`, and prepares input files for
the :class:`blast` and :class:`msa` modules with :class:`its_io.writer`.
"""

import os
import re
import pandas
import string
import random
import logging

from Bio import SeqIO
from os import path

from ab12phylo import filter


class reader:
    """
    The File Reader for :mod:`ABI` trace files in :code:`.ab1` format and :code:`.csv`
    tables of well-to-isolate coordinates. Generally reads in all :code:`.ab1` files in
    the provided directory, but accepts whitelist files for subsetting on sample IDs as
    well as files. Does some trimming, quality control and renaming using :class:`filter`.
    Sequences are saved in the per-gene dictionary :code:`self.seqdata` and source file
    names in :code:`self.metadata`. Reference sequences are renamed to avoid issues with
    :mod:`Gblocks` or the :class:`msa` tools. Sets :code:`args.genes` to the found genes
    if it was :code:`None`, favouring :code:`ITS1F`.

    :param args: :class:`argparse.Namespace` containing CLI and config information
    :returns None: object attributes are directly accessed
    """

    def __init__(self, args):
        """
        Reads in ABI trace files and .csv tables of well-to-isolate-coordinates
        using the arguments supplied to the package.

        :param args: argparse.Namespace containing CLI and config information
        :returns None: object attributes are directly accessed
        """
        self.log = logging.getLogger(__name__)
        self.args = args
        self.csvs, self.seqdata, self.metadata = dict(), dict(), dict()

        self._get_ids()
        self._get_seqs()

        # for empty genes list, use the ones that were read in
        if self.args.genes is None:
            self.args.genes = list(self.seqdata.keys())

            # prefer ITS as guiding gene. no ref-reordering necessary because no in by_order case
            if 'ITS1F' in self.args.genes:
                self.args.genes.insert(0, self.args.genes.pop(args.genes.index('ITS1F')))

        self.log.info('--GENES-- %s' % '::'.join(self.args.genes))  # leave it alone!

        self._get_refs()
        self.seq_counts = {gene: len(records) for gene, records in self.seqdata.items()}

    def _get_ids(self):
        """
        Reads in well-to-isolate coordinates from all .csv files in :param:`csv_dir`
        as :class:`pandas.DataFrame` objects. Box numbers are extracted from filenames.
        """
        for root, dirs, files in os.walk(self.args.csv_dir):
            for file in files:
                if file.endswith('.csv'):
                    df = pandas.read_csv(path.join(root, file), header=None)
                    df.index = list(range(1, df.shape[0] + 1))
                    df.columns = list(string.ascii_uppercase[0:df.shape[1]])
                    box = re.sub(r'\D', '', str(file).split('_')[-1])
                    self.csvs[box] = df
        if len(self.csvs) == 0:
            self.log.error('No .csv files found.')
            exit('No .csv files found.')
        self.log.debug('read %d .csv files in: %s' % (len(self.csvs), self.args.csv_dir))
        return

    def _get_seqs(self):
        """
        Reads in all .ab1 trace files in or below :param:`args.abi_dir` as
        :class:`Bio.SeqRecord` objects and renames to sample IDs from :param:`csvs`.
        Keeps track of source files, trims sequences and marks low-quality regions using
        filter.py module. Skips seq if file parsing or coordinate extraction fails, or
        sequence has too low quality. Versionize if ID is duplicated.
        """

        # Prepare subset analysis via file or sample ID
        file_subsetting = sample_subsetting = False
        abi_whitelist = list()
        if self.args.abi_set is not None:
            file_subsetting = True
            abi_whitelist = open(self.args.abi_set, 'r').read().strip().split('\n')
            # make absolute paths
            abi_whitelist = [path.join(path.dirname(path.dirname(__file__)), abi.strip()[1:])
                             if abi != '' and abi[0] == '$' else abi.strip() for abi in abi_whitelist]
            self.log.debug('read ABI trace file whitelist with %d entries' % len(abi_whitelist))

        sample_whitelist = list()
        if self.args.sample_set is not None:
            sample_subsetting = True
            sample_whitelist = open(self.args.sample_set, 'r').read().strip().split('\n')
            self.log.info('read sample ID whitelist with %d entries' % len(sample_whitelist))

        count = 0
        pattern = re.compile(r'FN_pl([\d])_[\d]{1,2}_(.*)_([A-Z][\d]{2})_')

        with open(self.args.bad_seqs, 'w') as bad_seqs:
            bad_seqs.write('file\tid\tbox\tgene\tproblem\n')

            for root, dirs, files in os.walk(self.args.abi_dir):
                for file in files:
                    if file.endswith('.ab1'):

                        # skip files not listed if using a restricted set
                        if file_subsetting and path.join(root, file) not in abi_whitelist:
                            continue

                        count += 1

                        # read ABI file
                        try:
                            record = SeqIO.read(path.join(root, file), 'abi')
                        except UnicodeDecodeError:
                            bad_seqs.write('%s\t \t \t \tSeqIO error\n' % file)
                            self.log.error('SeqIO error: %s' % file)
                            continue

                        # parse box, gene and coordinates from name
                        # NOTE: ID of a SeqRecord sometimes indicates wrong box. -> use name.
                        try:
                            m = pattern.search(record.name)
                            box, gene, coords = m.groups()
                            (row, col) = (int(coords[1:]), coords[0])
                            try:
                                # swap out well coordinates for isolate numbers
                                record.id = self.csvs[box].loc[row, col]
                            except KeyError:
                                bad_seqs.write('%s\t \t%s\t%s\tbox not found\n' % (file, box, gene))
                                self.log.warning('SeqIO box %s not found %s' % (box, file))
                                continue

                            # filter with whitelist
                            if sample_subsetting and record.id not in sample_whitelist:
                                continue

                            # skip seqs from wrong gene
                            if self.args.genes is not None and gene not in self.args.genes:
                                bad_seqs.write('%s\t%s\t%s\t%s\twrong gene\n' % (file, record.id, box, gene))
                                self.log.info('SeqIO wrong gene %s' % file)
                                continue

                        except AttributeError:
                            bad_seqs.write('%s\t \t \t \tunexpected filename\n' % file)
                            self.log.error('SeqIO unexpected name of ABI record %s' % file)
                            continue

                        # quality check
                        try:
                            record = filter.trim_ends(record, self.args.min_phred, self.args.end_ratio)
                            record = filter.mark_bad_stretches(record, self.args.min_phred, self.args.bad_stretch)
                        except ValueError:
                            bad_seqs.write('%s\t%s\t%s\t%s\tlow quality\n' % (file, record.id, box, gene))
                            self.log.warning('SeqIO low quality: %s' % file)
                            continue

                        # ensure gene dict is present
                        if gene not in self.seqdata:
                            self.seqdata[gene] = dict()
                            self.metadata[gene] = dict()

                        elif record.id in self.seqdata[gene]:
                            old_id = record.id
                            # add suffix to duplicate IDs
                            record = filter.new_version(record, self.seqdata[gene].keys())
                            bad_seqs.write('%s\t%s\t%s\t%s\tnew id, %s already present\n'
                                           % (file, record.id, box, gene, old_id))
                            self.log.info('new versioned id: %s from %s' % (record.id, file))

                        # save SeqRecord by position
                        self.seqdata[gene][record.id] = record

                        self.metadata[gene][record.id] = {'file': path.join(root, file), 'box': box}

        if len(self.seqdata) == 0:
            self.log.error('No .ab1 ABI trace files found.')
            exit('I/O: No .ab1 ABI trace files found.')
        self.log.debug('read %d .ab1 files in: %s' % (count, self.args.abi_dir))
        self.log.debug('read %d sequences of acceptable quality' %
                       sum(len(gene_records) for gene_records in self.seqdata.values()))
        return

    def _get_refs(self):
        """
        Reads in reference sequences from :code:`.fasta` reference file to guiding gene in
        the per-gene-per-box dictionary of seqs. IDs are replaced by short, random ids and
        original descriptions are saved in the :code:`metadata` dictionary.
        """
        if self.args.ref is None:
            return

        # use random but reproducible prefixes. different seed from raxml -> less confusion
        random.seed(self.args.seed + 2)

        # matching reference files by name to genes
        if not self.args.by_order:
            message = 'Trying to match genes and reference file names. If allocation is wrong, ' \
                      'try setting --ref and --genes manually, with ITS1F possibly first argument.'
            # convert all genes to lowercase
            genes = [gene.lower() for gene in self.args.genes]
            # extract filenames from refs
            extracts = [ref[:-6].split('/')[-1].lower() for ref in self.args.ref]

            order = list()
            for filename in extracts:
                try:
                    order.append(genes.index(filename))
                except ValueError:
                    order.append(-1)
        else:
            message = 'Assigning references by order. If mismatched, re-order --ref arguments accordingly.'
            order = list(range(len(self.args.ref)))

        self._warn_order(message, order)

        # dict of reference source organism -> random key
        lookup = dict()

        # read in refs
        for ref_file, gene_pos in zip(self.args.ref, order):
            gene = self.args.genes[gene_pos]
            count = 0

            for record in SeqIO.parse(ref_file, 'fasta'):
                count += 1
                # parse species and possibly strain
                strain = re.split(r'[\s_]', record.description.strip().split(',')[0])
                accession = strain.pop(0)
                try:
                    strain = ' '.join(strain[:strain.index('strain') + 3])
                except ValueError:
                    strain = ' '.join(strain[:3])

                # retrieve or generate random key
                if strain in lookup:
                    _id = lookup[strain]
                else:
                    # swap original ids for random short ones that no tool takes offense at.
                    _id = 'REF_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
                    lookup[strain] = _id

                # save original id+description
                self.metadata[gene][_id] = {'file': ref_file,
                                            'accession': accession,
                                            'reference_species': strain}
                record.id = _id
                record.description = ''  # MARK do not delete deletion
                # save ref record
                self.seqdata[gene][_id] = record

            self.log.debug('Read %d reference sequences for %s from: %s'
                           % (count, gene, ref_file))
        return

    def _warn_order(self, message, order):
        """
        Builds a table of gene <- ref allocations. If there are indications of errors, logs as warning.
        :param order: positions of genes that references were matched to
        """
        table = 'allocation table\n\t<gene>\t<-\t<reference>'
        for j in range(len(self.args.genes)):
            table += '\n\t%s\t<-\t' % self.args.genes[j]
            try:
                # the reference where the order is j
                table += self.args.ref[order.index(j)]
            except ValueError:
                table += '___no ref___'

        # log according to error indications
        if self.args.by_order and len(order) > 1 or -1 in order or len(order) != len(self.args.genes):
            self.log.warning(message)
            self.log.warning(table)
        else:
            self.log.info(message)
            self.log.info(table)

        # pass over order anyway
        for i in range(len(order)):
            if order[i] == -1:
                self.log.warning('Reference file could not be matched to a gene: %s' % self.args.ref[i])

        # pass over genes
        for j in range(len(self.args.genes)):
            if j not in order:
                self.log.warning('Gene %s has no reference.' % self.args.genes[j])
        return


class writer:
    """
    The :mod:`FASTA` File Writer. For each gene, :code:`<results>/<gene>/<gene>.fasta`
    will contain all sequences for the :class:`blast` and :class:`msa` modules.

    :param args: :class:`argparse.Namespace` containing CLI and config information
    :param _reader: :class:`its_io.reader` with seq and meta data
    :returns None:
    """

    def __init__(self, args, _reader):
        self.log = logging.getLogger(__name__)

        for gene in _reader.seqdata:
            # make a subdir for <gene> data
            os.makedirs(path.join(args.dir, gene), exist_ok=True)

            # write FASTA
            with open(path.join(args.dir, gene, gene + '.fasta'), 'w') as fasta:
                SeqIO.write(_reader.seqdata[gene].values(), fasta, 'fasta')
                self.log.debug('wrote %s seqdata %s' % (gene, fasta.name))

        # make a pandas DataFrame from metadata
        df = pandas.concat({gene: pandas.DataFrame.from_dict(records, orient='index')
                            for gene, records in _reader.metadata.items()})
        # write to TSV
        df.to_csv(args.tsv, sep='\t', na_rep='', header=True, index=True)
        self.log.debug('wrote metadata to %s' % args.tsv)
