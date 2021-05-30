#!/usr/bin/env python3
# 2021 Leo Kaindl

import logging
import os
import sys

LOG = logging.getLogger(__name__)

try:
    import gi

    gi.require_version('Gtk', '3.0')
except ModuleNotFoundError as ex:
    LOG.exception(ex)
    from subprocess import run, PIPE
    import shlex
    import shutil
    import getpass

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
            proc = run(shlex.split(arg))
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
            proc = run(shlex.split(apt_call))
            if proc.returncode != 0:
                exit(1)
        except Exception as ex:
            LOG.exception(ex)

from ab12phylo import repo
from ab12phylo.gtk_blast import blast_page
from ab12phylo.gtk_gbl import gbl_page
from ab12phylo.gtk_io import io_page
from ab12phylo.gtk_ml import ml_page
from ab12phylo.gtk_msa import msa_page
from ab12phylo.gtk_qal import qal_page
from ab12phylo.gtk_rgx import rgx_page
from ab12phylo.gtk_tree import tree_page


# set the icon theme
# Gtk.Settings.get_default().set_property('gtk-icon-theme-name', 'Papirus-Dark-Maia')
# Gtk.Settings.get_default().set_property('gtk-theme-name', 'Matcha-dark-sea')


class ab12phylo_app(io_page, rgx_page, qal_page, msa_page,
                    gbl_page, blast_page, ml_page, tree_page):

    def __init__(self):
        super().__init__()
        self.supers = [io_page, rgx_page, qal_page, msa_page,
                       gbl_page, blast_page, ml_page, tree_page]
        self.re_runs = {2: self.start_trim,  # 1: self.start_read,
                        4: self.start_gbl, 7: self.start_phy}

    def load(self, path, *args):
        super().load(path)
        for module in self.supers:
            try:
                module.reload_ui_state(self)
            except AttributeError as e:
                LOG.debug(e)

    def reload_ui_state(self, page=False, *args):
        if not page:
            return
        try:
            self.supers[page].reload_ui_state(self)
        except AttributeError as e:
            LOG.debug(e)

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        # tell the MSA pre-set about it
        if 'aligner' in self.iface:
            self.refresh_paths()

    def test(self, *args):
        super().test(self, args)
        self.supers[1].refresh(self)

    def refresh(self, page=-1, *args):
        """
        Call the refresh function of the current page and hide or show
        the appropriate buttons.
        """
        iface = self.iface
        if page == -1 or type(page) != int:
            page = iface.notebook.get_current_page()
        self.supers[page].refresh(self)
        # hide or show these two actions depending on applicability
        iface.refresh.props.visible = bool(page in self.re_runs)
        iface.reset.props.visible = page in {1, 2, 4, 7}
        iface.gene_roll.props.visible = page == 2
        iface.next.props.visible = page != 7
        iface.back.props.visible = page != 0
        iface.open_test.props.visible = page == 0
        iface.helper.set_markup(repo.help.get(page, ''))

    def re_run(self, *args):
        """
        Depending on the currently visible page, re-run the matching background task.
        Handles the Refresh button.
        """
        page = self.iface.notebook.get_current_page()
        if page in self.re_runs:
            self.re_runs[page]()
        else:
            raise RuntimeWarning(f'page {page} has no re-run')

    def reset(self, *args):
        """
        Depending on the currently visible page (either rgx, trim, gbl or tree), reset the
        regex table or the tree. Handles the Reset button.
        """
        page = self.iface.notebook.get_current_page()
        if page == 1:
            self.set_changed(1)
            self.reset_columns(do_parse=True)
        elif page == 2:
            self.reset_x_scale()
            self.data.qal.ignore_ids = {g: set() for g in self.data.genes}
            self.start_read(run_after=[self.start_trim])
            # self.start_trim()
        elif page == 4:
            self.reset_x_scale()
            self.undrop_seqs()
        elif page == 7:
            self.reset_tree()

    def proceed(self, widget=None, page=None):
        """
        The function connected to the _Next button. For pages with a background thread,
        this will start it and instruct it to re-run this function afterwards; to make the
        application proceed only upon thread completion.
        :param widget: optional, for use as callback
        :param gui:
        :param page:
        :return:
        """
        data = self.data
        iface = self.iface
        page = iface.notebook.get_current_page() if not page else page

        if self.get_errors(page):
            self.show_notification('There are still errors on the page!')
            return

        if self.get_changed(page):
            if page == 0:
                self.supers[page].refresh(self)
                self.reset_columns(do_parse=True)
            elif page == 1:
                self.supers[page].refresh(self)
                if 1 < sum(data.rx_fired) < 5:
                    self.show_notification('Make sure all columns have been parsed.')
                    return
                self.start_read(run_after=[self.init_gene_roll, self.proceed])
                return  # leave this alone
            elif page == 2:
                self.trim_all(run_after=[self.proceed])
                return
            elif page == 3:
                self.start_align(None, run_after=[self.proceed])
                return
            elif page == 4:
                self.start_gbl()
            elif page == 5:
                iface.tempspace.df.to_csv(
                    self.wd / repo.PATHS.tsv, sep='\t', na_rep='', header=True, index=True)
                self.set_changed(5, False)
            elif page == 6:
                self.start_phy()
            self.set_changed(page, False)

        # then proceed
        iface.notebook.next_page()
        data.page = iface.notebook.get_current_page()
        LOG.debug('proceeded to page %d' % data.page)

    def tree_modify(self, action, sel_idx=None):
        """
        Fetch the names of the selected nodes.
        :param action:
        :param sel_idx:
        :return: {action: [matching names]} that will be added to the existing modifications
        """
        if not sel_idx:
            sel_idx = [tp.get_indices()[0] for tp in self.iface.tree_sel.get_selected_rows()[1]]
        if not sel_idx:
            self.show_notification('Cannot %s, no selection' % action.get_name(), secs=2)
            return
        if not self.iface.tree:
            self.show_notification('No tree in memory, re-drawing first', secs=2)
            self.start_phy(run_after=(self.tree_modify, (action, sel_idx)))
            return

        # return a list of tip labels
        idsp = [(_id, sp) for i, (_id, sp) in enumerate(
            self.data.tree_anno_model.get_column((0, 1))) if i in sel_idx]

        names = list()
        for _id, sp in idsp:
            if sp.startswith('REF_'):
                names.append(sp.split(':')[0])
            else:
                names.append(_id.replace(', ', '~'))
        LOG.debug('re-plot and %s' % action.get_name())
        self.start_phy({action.get_name(): names})

    def delete_and_ignore_rows(self, *args):
        if super().delete_and_ignore_rows(*args):
            page = args[2]
            if page == 2:
                self.iface.view_qal.grab_focus()
                self.start_trim()
            elif page == 4:
                self.iface.view_gbl.grab_focus()
                self.drop_seqs()
        return True


def main():
    app = ab12phylo_app()
    exit_status = app.run(sys.argv)
    sys.exit(exit_status)


if __name__ == '__main__':
    main()
