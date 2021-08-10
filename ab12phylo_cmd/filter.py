# 2021 Leo Kaindl

import configparser
import copy
import os
import re
import shutil
import stat
import sys
import tarfile
import zipfile
from pathlib import Path
from subprocess import run, check_output, PIPE

import requests
from Bio.Seq import MutableSeq
from bs4 import BeautifulSoup

regex = re.compile('\\.[\\d]+$')


def trim_ends(seqrecord, min_phred, end_ratio, trim_preview=False):
    """
    Trims a supplied SeqRecord with given minimal quality score and ratio of good bases in end of given length.
    """
    if not seqrecord.letter_annotations:
        raise AttributeError('no quality')
    phreds = seqrecord.letter_annotations['phred_quality']
    if set(phreds) == {0} and min_phred > 0:
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
    start = max(0, start)

    # trim right
    _ok = 0
    end = len(phreds)
    while _ok < end_ratio[0] and end >= 0:
        Ns = len([i for i in phreds[end - end_ratio[1]: end] if i < min_phred])
        end += -Ns + end_ratio[1] - end_ratio[0]
        _ok = end_ratio[1] - Ns

    if start >= end:
        raise ValueError('low quality')

    if trim_preview:
        # make some copies
        qal = copy.copy(phreds)
        rec = copy.copy(seqrecord)
        rec.letter_annotations = copy.copy(rec.letter_annotations)
        rec.letter_annotations['phred_quality'] = qal

        if type(rec.seq) != MutableSeq:
            rec.seq = rec.seq.tomutable()

        if start > 0:
            rec.seq[0:start] = '-' * start
            qal[0:start] = [min_phred + 1] * start

        assert end > 0
        rec.seq[end:] = '-' * (len(rec) - end)
        qal[end:] = [min_phred + 1] * (len(rec) - end)
        return rec
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

    if not seqrecord.letter_annotations:
        raise AttributeError('no quality')
    phreds = seqrecord.letter_annotations['phred_quality']
    if type(seqrecord.seq) != MutableSeq:
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


def new_id(rid, keys):
    while rid in keys:
        match = regex.search(rid)
        if match:
            rid = rid[:match.start()] + '.' \
                  + str(int(match.group()[1:]) + 1)
        else:
            rid += '.1'
    return rid


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


def chmod_x(p):
    """Make a file executable. Accepts either a string or a pathlib.Path"""
    if type(p) == str:
        p = Path(p).resolve()
    try:
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    except FileNotFoundError:
        pass


def fetch_non_python_tools(suffix, cfg, other_cfg, save_dir, log):
    """
    :param suffix: Indicator it this is ab12phylo -> '' or ab12phylo-cmd -> '-cmd'
    :param cfg: Path to the 'conf.cfg' file
    :param other_cfg: Path to the 'conf.cfg' file for the other AB12PHYLO version
    :param save_dir: Installation target directory. Actually always ./tools, where . is where this file is
    :param log: logging.Logger instance
    :return:
    """
    paths = dict()
    all_tools = ['blastn', 'raxml-ng', 'iqtree2']
    skip = False

    def find(tool, fix_linux=False):
        try:
            exe_suffix = '.exe' if sys.executable.endswith('.exe') and not fix_linux else ''
            exe = next(save_dir.rglob(f'{tool}{exe_suffix}'))
            # Make the file executable
            exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
            log.info(f'Found {tool} at {exe}')
            if not fix_linux:
                paths[tool] = exe
            else:
                paths[tool + '-linux'] = exe
        except StopIteration:
            pass

    def prompt(tool, source, alt_text=False):
        answer = ''
        while answer not in {'y', 'yes', 'n', 'no'}:
            if not alt_text:
                answer = input(f'Download {tool} from {source}? [y/n]').lower().strip()
            else:
                answer = input('Import paths and skip re-download? [y/n]').lower().strip()
        return answer in {'y', 'yes'}

    if other_cfg.is_file():
        other_paths = dict()
        # fetch the paths from the config
        cfg_parser = configparser.ConfigParser()
        cfg_parser.read(other_cfg)
        if 'Paths' in cfg_parser:
            other_paths.update(dict(cfg_parser['Paths']))
            if all(tool in other_paths for tool in all_tools):
                # if all three paths are there, offer to skip the re-download
                log.info(f'Found config for other AB12PHYLO version at {other_cfg}')
                if prompt('', '', True):
                    skip = True  # The user chose to not download new tools
                    paths.update(other_paths)

    for tool in all_tools:
        if skip:
            continue

        # Look for a local installation and check if it is recent enough
        exe = shutil.which(tool)
        if exe is not None:
            log.info(f'Found {tool} at {exe}')
            if tool == 'blastn':
                output = check_output(exe + ' -version', shell=True).decode('utf-8')
                version = [int(i) for i in output.split(',')[0].split(' ')[-1].split('.')]
                if version[0] == 2 and version[1] >= 9 or version[0] > 2:
                    paths[tool] = exe
                    continue
                else:
                    log.warning('installed BLAST+ is outdated')

            elif tool == 'raxml-ng':
                # check if the installed version is old-ish / the packaged version is newer
                proc = run(stdout=PIPE, stderr=PIPE, shell=True, args=f'{exe} -v')
                version = [int(i) for i in re.search(
                    'RAxML-NG v[\.\s]+?([0-9\.]+)', proc.stdout.decode()).groups()[0].split('.')]
                if version[0] == 0 or version[0] == 1 and version[1] == 0 and version[2] <= 1:
                    log.warning('installed RAxML-NG is outdated')
                else:
                    paths[tool] = exe
                    continue

        # No or outdated local installation
        # At this point, always doing a fresh download -> enables updating.
        try:
            if tool == 'blastn':
                log.info(f'{tool} not installed (not on the $PATH).')
                if not prompt('BLAST+', 'the NCBI'):
                    continue  # The user chose to not download a new BLAST+
                url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST'
                log.debug(f'Fetching {url} ... ')
                r = requests.get(url, timeout=12)
                html = BeautifulSoup(r.content, 'html.parser')
                links = [tag['href'] for tag in html.find_all('a') if 'href' in tag.attrs]

                archive_suffix = '-x64-%s.tar.gz' % {'linux': 'linux', 'win32': 'win64',
                                                     'darwin': 'macosx', 'cygwin': 'win64'}[sys.platform]

                right_one = [l for l in links if l.endswith(archive_suffix)][0]
                zf = save_dir / right_one

                log.debug(f'Downloading {right_one} ... ')
                r = requests.get(f'{url}/{right_one}', timeout=12)
                with open(zf, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=128):
                        fd.write(chunk)

                log.debug(f'Extracting BLAST+ ... ')
                with tarfile.open(zf) as zo:
                    zo.extractall(zf.parent)
                zf.unlink()
                find(tool)

            elif tool in ['raxml-ng', 'iqtree2']:
                if not prompt({'raxml-ng': 'RAxML-NG',
                               'iqtree2': 'IQ-Tree'}[tool], 'GitHub'):
                    continue

                platform = sys.platform
                url = 'https://api.github.com/repos/%s/%s/releases/latest'
                if tool == 'raxml-ng':
                    url %= 'amkozlov', 'raxml-ng'
                else:
                    url %= 'iqtree', 'iqtree2'
                log.debug(f'Fetching {url} ... ')
                js = requests.get(url, timeout=12).json()
                for asset in js['assets']:
                    name, path = asset['name'], asset['browser_download_url']
                    # find the right release, but make sure to always download a linux version
                    if not (tool == 'raxml-ng' and ('linux_x86_64.zip' in name
                                                    or platform == 'darwin'
                                                    and 'macos_x86_64.zip' in name)
                            or tool.startswith('iqtree2') and ('Linux.tar.gz' in name
                                                               or platform == 'win32'
                                                               and 'Windows.zip' in name
                                                               or platform == 'darwin'
                                                               and 'MacOSX.zip' in name)):
                        continue
                    # download
                    log.debug(f'Downloading {name} ... ')
                    r = requests.get(path, stream=True)
                    zf = save_dir / name
                    with open(zf, 'wb') as file:
                        for chunk in r.iter_content(chunk_size=128):
                            file.write(chunk)

                    # extract
                    log.debug(f'Extracting {zf} ... ')
                    zs = zf.with_suffix('')
                    if zf.suffix == '.zip':
                        with zipfile.ZipFile(zf, 'r') as zo:
                            # prevent ridiculous double level
                            if all(n.startswith(zs.name) for n in zo.namelist()):
                                zo.extractall(zs.parent)
                            else:
                                zo.extractall(zs)
                    elif zf.suffixes[-2:] == ['.tar', '.gz']:
                        zs = zs.with_suffix('')  # yes, again
                        with tarfile.open(zf) as zo:
                            # prevent ridiculous double level
                            if all(n.startswith(zs.name) for n in zo.getnames()):
                                zo.extractall(zs.parent)
                            else:
                                zo.extractall(zs)
                    zf.unlink()
                    if platform == 'linux' or 'Linux.tar.gz' not in name:
                        # list the tool for non-linux platforms under its name
                        find(tool)
                    else:
                        # list the linux version with the '-linux' suffix
                        find(tool, fix_linux=True)

        except requests.ConnectionError as ex:
            log.error('Couldn\'t download %s. Is the system offline?' % tool)
            log.error(ex)
        except Exception as ex:
            log.error('Downloading %s failed.' % tool)
            log.error(ex)

    # when nothing was downloaded, check if the tools are there and save paths
    for tool in all_tools:
        if tool not in paths:
            find(tool)

    with open(cfg, 'w') as fh:
        fh.write(f'''
# This file indicates that ab12phylo{suffix} has been run before
# and stores paths to some non-python external tools.
# Invoking ab12phylo{suffix} with the `--initialize` flag or 
# deleting the file will trigger a re-download of these 
# tools the next time ab12phylo{suffix} is run.\n''')
        fh.write(f'\n[Paths]')
        for k, v in paths.items():
            fh.write(f'\n{k}={v}')
