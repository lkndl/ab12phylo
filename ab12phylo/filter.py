from Bio.Seq import MutableSeq


def trim_ends(seqrecord, min_phred, end_ratio, keep=False):
    """
    Trims a supplied SeqRecord with given minimal quality score and ratio of good bases in end of given length.
    """
    if not seqrecord.letter_annotations:
        raise AttributeError('no quality')
    phreds = seqrecord.letter_annotations['phred_quality']
    if set(phreds) == {0}:
        raise AttributeError('no quality')

    # trim left
    start = _ok = 0
    while _ok < end_ratio[0] and start < len(phreds):
        Ns = len([i for i in phreds[start: start + end_ratio[1]] if i < min_phred])
        # assume all good bases are at right edge of counting window
        start += Ns
        # and allow for some leading non-confident bases
        start -= (end_ratio[1] - end_ratio[0])
        _ok = end_ratio[1] - Ns

    # trim right
    _ok = 0
    end = len(phreds)
    while _ok < end_ratio[0] and end >= 0:
        Ns = len([i for i in phreds[end - end_ratio[1]: end] if i < min_phred])
        end += -Ns + end_ratio[1] - end_ratio[0]
        _ok = end_ratio[1] - Ns

    if start >= end:
        raise ValueError('low quality')

    if keep:
        seqrecord.seq = seqrecord.seq.tomutable()
        seqrecord.seq[0:start] = '-' * start
        seqrecord.seq[end:] = '-' * (len(seqrecord) - end)
        return seqrecord
    else:
        return seqrecord[start:end]


def mark_bad_stretches(seqrecord, min_phred, _len):
    """
    Replaces stretches of bases that were sequenced with below-minimal phred quality score in a
    :class:`Bio.SeqRecord` object by a sequence of undetermined (:code:`N`) nucleotides of equal length.
    *Bad stretches* shorter than :code:`_len` are not changed.

    :param seqrecord: The `SeqRecord` object that will be filtered.
    :param min_phred: The minimal score below with replacement with ``N`` might occur
    :param _len: The minimal length of bad bases that will be replaced.
    :return: The input ``SeqRecord``, with bad stretches replaced.
    """

    phreds = seqrecord.letter_annotations['phred_quality']
    seqrecord.seq = seqrecord.seq.tomutable()
    start = 0

    while True:
        # look for bad base
        if phreds[start] < min_phred:
            end = start + 1
            while end < len(phreds) and phreds[end] < min_phred:
                # elongate bad stretch
                end += 1
            if end - start >= _len:
                # replace long bad stretch
                seqrecord.seq[start: end] = 'N' * (end - start)
            start += 1
        else:
            # look elsewhere
            start += 1
        if start >= len(seqrecord):
            break

    return seqrecord


def new_version(seqrecord, keys):
    """Updates duplicate isolate coordinates with increasing version number."""
    while seqrecord.id in keys:
        pos = seqrecord.id.find('.')
        if pos == -1:
            seqrecord.id += '.1'
        else:
            version = int(seqrecord.id[pos + 1:])
            seqrecord.id = seqrecord.id[:pos + 1] + str(version + 1)
    return seqrecord


def mark_bad_bases(seqrecord, min_phred):
    """
    unused: mark *all* non-confident bases as N
    """
    phreds = seqrecord.letter_annotations['phred_quality']
    seq = str()
    for pos in range(len(phreds)):
        if phreds[pos] < min_phred:
            seq += 'N'
        else:
            seq += seqrecord.seq[pos]
    seqrecord.seq = MutableSeq(seq)
