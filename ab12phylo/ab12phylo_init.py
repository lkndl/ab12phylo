#!/usr/bin/env python3
# 2021 Leo Kaindl

import logging
import os
import sys

LOG = logging.getLogger(__name__)


def check_pygobject():
    """
    Check if GTK3 is functional, meaning if 'gi' can be imported and is recent enough.
    If this python is managed by conda, suggest and run a conda install of pygobject, gtk3,
    and the missing icons. Alternatively and if the package manager apt is available,
    suggest a shot-in-the-dark apt install command.
    :return:
    """
    try:
        import gi

        gi.require_version('Gtk', '3.0')
    except (ImportError, ModuleNotFoundError, ValueError) as ex:
        LOG.exception(ex)

        from subprocess import run, PIPE
        import shutil
        import getpass
        import os
        os.system('color')  # enable color console on windows

        yes = {'y', 'yes'}
        yes_or_no = yes.union({'n', 'no'})

        # anaconda installation
        if os.path.exists(os.path.join(sys.prefix, "conda-meta")):
            run_0 = f'conda activate {sys.prefix.split("/")[-1]}'
            run_1 = 'conda install -y -c conda-forge pygobject gtk3 '
            run_2 = 'conda install -y -c conda-forge adwaita-icon-theme hicolor-icon-theme '
            print(f'\nPyGObject was not found on your python at '
                  f'{sys.prefix}. It looks like you are using anaconda, '
                  f'so the missing dependencies can be installed via:\n'
                  f'\033[94m{run_0}\n{run_1}\n{run_2}\033[0m\n')

            if sys.prefix.endswith('base'):
                print(f'\033[91mInstalling to the \033[0mbase\033[91m environment '
                      f'is likely to cause package conflicts!\033[0m\nYou can find '
                      f'tested environments for Linux https://raw.githubusercontent'
                      f'.com/lkndl/ab12phylo/main/recipe/ab1.yaml and Windows https'
                      f'://raw.githubusercontent.com/lkndl/ab12phylo/main/recipe/wi'
                      f'n.yaml in the AB12PHYLO repo. Install one with:\n'
                      f'\033[94mconda env create -f ab1.yaml\n'
                      f'conda activate ab1\033[0m\n')

            answer = ''
            while answer not in yes_or_no:
                answer = input(f'Install now? [y/n]').lower().strip()
            if answer in {'n', 'no'}:
                exit(1)

            print(f'\033[94m{sys.prefix}\033[0m')
            for arg in [run_1, run_2]:
                print(f'\033[94m{arg}\033[0m')
                proc = run(arg, shell=True)
                if proc.returncode != 0:
                    exit(1)
            try:
                import gi

                gi.require_version('Gtk', '3.0')
            except Exception as ex:
                LOG.exception(ex)
        else:
            # non-anaconda case
            print(f'\nPyGObject was not found on your python at '
                  f'{sys.prefix}.')
            if shutil.which('apt') is None:
                print('On this system, please install AB12PHYLO to a python3 '
                      'environment managed by the anaconda package manager.')
                exit(1)

            apt_call = 'sudo apt install python3-gi python3-gi-cairo gir1.2-gtk-3.0'
            print(f'It looks like you are using a Linux distro '
                  f'that uses the Debian package manager APT. '
                  f'GTK3 is often already installed on such a system, '
                  f'but you can try installing the missing packages '
                  f'with\n\033[94m{apt_call}\033[0m\nbut this might break '
                  f'your system!\nA safer method would be using the '
                  f'anaconda package manager and setting up a tested '
                  f'environment: https://raw.githubusercontent.com/lkndl/'
                  f'ab12phylo/main/recipe/ab1.yaml. Create it via:\n'
                  f'\033[94mconda env create -f ab1.yaml\n'
                  f'conda activate ab1\033[0m\n')
            answer = ''
            while answer not in yes_or_no:
                answer = input(f'Try apt installation now? \033[91mUse with '
                               f'caution!\033[0m [y/n]').lower().strip()
            if answer in {'n', 'no'}:
                exit(1)

            print(f'\033[94m{apt_call}\033[0m')
            try:
                proc = run(apt_call, shell=True)
                if proc.returncode != 0:
                    exit(1)
            except Exception as ex:
                LOG.exception(ex)


def initialize():
    """
    For Linux systems, create a desktop file.
    Try downloading the test data set.
    Trigger searching for BLAST+, RAxML-NG and iqtree2 installations via shutil.
    If a tool is not found on the $PATH or outdated, run prompts and installations.
    :return:
    """
    import requests
    from pathlib import Path

    from ab12phylo import repo
    from ab12phylo.__init__ import __version__
    from ab12phylo_cmd.filter import fetch_non_python_tools

    LOG.info('Downloading test data ...')
    r = requests.get('https://github.com/lkndl/ab12phylo/wiki/test_data.zip',
                     stream=True, timeout=20)
    zf = repo.BASE_DIR / 'ab12phylo' / 'test_data.zip'
    with open(zf, 'wb') as file:
        for chunk in r.iter_content(chunk_size=128):
            file.write(chunk)
    # leave as a ZIP!

    if sys.platform in {'linux', 'darwin'}:
        LOG.info('Try creating a desktop entry ...')
        try:
            for dsk in [Path('~/.local/share/applications'), Path('/usr/local/share/applications')]:
                dsk = dsk.expanduser()
                if not dsk.is_dir():
                    continue
                with open(dsk / 'ab12phylo.desktop', 'w') as file:
                    file.write('''[Desktop Entry]
Name=AB12PHYLO
Version={version}
Comment=Integrated pipeline for ML phylogenetic inference from ABI trace and FASTA data
Exec={py} {script}
Icon={icon}
Terminal=false
Type=Application
Categories=GTK;GNOME;Utility;Application;
StartupNotify=true'''.format(version=__version__, py=sys.executable,
                             script=repo.BASE_DIR / 'ab12phylo' / 'ab12phylo_app.py',
                             icon=repo.PATHS.icon_path.with_suffix('.svg')))
                break

        except Exception as ex:
            LOG.info(ex)

    fetch_non_python_tools('', repo.BASE_DIR / 'ab12phylo' / 'conf.cfg',
                           repo.PATHS.cmd_config, repo.TOOLS, LOG)


def main():
    LOG.info('Initializing ab12phylo:')
    check_pygobject()
    initialize()


if __name__ == '__main__':
    main()
