# 2021 Leo Kaindl

import logging
import os
import random
import re
import shlex
import shutil
import stat
import threading
from pathlib import Path
from subprocess import run, Popen, PIPE
from time import sleep, time
from zipfile import ZipFile, ZIP_DEFLATED

import gi

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo import repo
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 6


class ml_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.evo_model.set_model(data.evo_models)
        iface.evo_model.get_child().connect(
            'focus_out_event', self._change_evo_model, iface.evo_modify, data.ml)
        iface.evo_block = iface.evo_model.connect_after(
            'changed', self._load_model_file, iface.evo_modify)
        iface.evo_model.handler_block(iface.evo_block)

        iface.raxml_run.connect('clicked', self.start_ML, 'raxml')
        iface.raxml_export.connect('clicked', self.start_ML, 'raxml_export')
        iface.raxml_import.connect('clicked', self.import_tree)

        for wi in [iface.bootstraps, iface.rand, iface.pars, iface.raxml_seed]:
            wi.connect('key-press-event', self.edit_numerical_entry_up_down)
            wi.connect('changed', self.edit_numerical_entry)

        iface.raxml_seen = False

    def import_tree(self, widget):
        """Import a tree file or two"""
        dialog = Gtk.FileChooserDialog(title='import tree(s)',
                                       parent=None, select_multiple=True,
                                       action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            paths = [Path(p).resolve() for p in dialog.get_filenames()]
            if len(paths) > 2:
                self.show_message_dialog('Please select at most two tree files.')
            else:
                errors = list()
                for path in paths:
                    if 'FBP' in path.name.upper():
                        shutil.copy(path, self.wd / repo.PATHS.fbp)
                    elif 'TBE' in path.name.upper():
                        shutil.copy(path, self.wd / repo.PATHS.tbe)
                    else:
                        errors.append(path.name)
                if errors:
                    self.show_message_dialog('Not immediately recognized as either '
                                             'tree_FBP.nwk or tree_TBE.nwk. You can also copy '
                                             'it/them to %s manually.' % self.wd, items=errors)
                else:
                    self.show_notification('imported trees:',
                                           [path.name for path in paths], 2)

        dialog.destroy()
        self.set_changed(PAGE, False)

    def _change_evo_model(self, entry, event_focus, evo_modify, ml):
        combo = entry.get_parent().get_parent()
        if combo.get_active_iter():
            ml.evo_model = entry.get_text()
            evo_modify.props.sensitive = True
        else:
            tx = entry.get_text()
            combo.get_model().append([tx, None])
            ml.evo_model = tx  # save in project dataset
            LOG.debug('entered custom evo model %s' % tx)

    def _load_model_file(self, combo, evo_modify):
        if combo.get_active_iter() and combo.get_active_id() == 'from file':
            dialog = Gtk.FileChooserDialog(title='select file with partition table',
                                           parent=None, select_multiple=False,
                                           action=Gtk.FileChooserAction.OPEN)
            dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                               Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
            response = dialog.run()
            if response == Gtk.ResponseType.OK:
                tx = Path(dialog.get_filename()).resolve()
                try:
                    Path.mkdir(self.wd / 'RAxML', exist_ok=True)
                    shutil.copy(tx, self.wd / 'RAxML' / 'user_model')
                except FileNotFoundError:
                    self.show_notification('File does not exist', secs=1)
                    dialog.destroy()
                    return
                except shutil.SameFileError:
                    pass
                evo_modify.props.sensitive = False
                combo.get_model().append([tx.name, str(tx)])
                combo.set_active_id(tx.name)
                LOG.debug('selected partitioned model file: %s -> RAxML/user_model' % str(tx))
            dialog.destroy()

    def _load_raxml_help(self):
        iface = self.iface
        binary = shutil.which('raxml-ng')
        if binary:
            # check if the packaged version is newer
            try:
                proc = run(stdout=PIPE, stderr=PIPE, shell=True, args='%s -v' % binary)
                version = [int(i) for i in re.search(
                    'RAxML-NG v[\.\s]+?([0-9\.]+)', proc.stdout.decode()).groups()[0].split('.')]
                # version = proc.stdout.decode().strip().split(
                #     '\n')[0].split('v.')[1].split('released')[0].strip()
                if version[0] == 0 or version[0] == 1 and version[1] == 0 and version[2] <= 1:
                    binary = None
                    LOG.info('prefer out-of-the-box RAxML-NG')
            except Exception as ex:
                LOG.error('RAxML version check failed')
                LOG.error(ex)
        if not binary:
            binary = str(repo.TOOLS / repo.PATHS.RAxML)
        iface.raxml_exe.set_filename(binary)
        res = run(stdout=PIPE, stderr=PIPE, shell=True,
                  args='%s --help' % binary)
        iface.ml_help.props.buffer.props.text = res.stdout.decode().lstrip()
        LOG.debug('got RAxML --help')

    def reload_ui_state(self):
        data = self.data
        iface = self.iface
        ml = data.ml
        for w_name in ['bootstraps', 'rand', 'pars']:
            iface.__getattribute__(w_name).set_text(str(ml.__getattribute__(w_name)))
        iface.raxml_seed.set_text(str(ml.raxml_seed) if 'raxml_seed' in ml else '')
        iface.evo_model.set_active_id(ml.evo_model)
        iface.evo_modify.set_text(ml.evo_modify)
        iface.raxml_shell.set_active(ml.raxml_shell)

    def refresh(self):
        """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
        LOG.debug('ML refresh')
        data = self.data
        iface = self.iface
        if not iface.raxml_seen:
            self._load_raxml_help()
            self.reload_ui_state()
            iface.evo_model.handler_unblock(iface.evo_block)
            iface.evo_model.set_active_id(data.ml.evo_model)
            iface.raxml_seen = True

        self.set_errors(PAGE, not any((self.wd / a).is_file()
                                      for a in [repo.PATHS.tbe, repo.PATHS.fbp]))

    def start_ML(self, widget, mode, run_after=None):
        """Set-up the ML inference thread"""
        data = self.data
        iface = self.iface
        ml = data.ml

        if not data.genes or iface.thread.is_alive():
            self.show_notification('Busy', secs=1)
            return

        if mode in ['raxml', 'raxml_export']:
            ml.raxml = iface.raxml_exe.get_filename()
            for w_name in ['evo_modify', 'bootstraps', 'rand', 'pars', 'raxml_seed']:
                ml.__setattr__(w_name, iface.__getattribute__(w_name).get_text())
                if w_name in ['bootstraps', 'rand', 'pars']:
                    ml.__setattr__(w_name, int(ml.__getattribute__(w_name)))
            ml.evo_model = data.evo_models[iface.evo_model.get_active()]
            if ml.evo_model[1]:
                ml.evo_modify = ''
                ml.evo_model = str(self.wd / 'RAxML' / 'user_model')
            else:
                ml.evo_model = ml.evo_model[0]

            ml.raxml_seed = random.randint(0, max(1000, ml.bootstraps)) \
                if ml.raxml_seed == '' else int(ml.raxml_seed)
            iface.raxml_seed.props.text = str(ml.raxml_seed)
            Path.mkdir(self.wd / 'RAxML', exist_ok=True)
            ml.raxml_shell = iface.raxml_shell.get_active()

            iface.run_after = run_after
            ml.prev = 0
            ml.key = False
            ml.stdout = list()
            ml.seen = {'ML': set(), 'BS': set()}
            ml.motifs = {'ML': 'ML tree search #', 'BS': 'Bootstrap tree #'}
            iface.k = ml.bootstraps + ml.rand + ml.pars + 3 if mode == 'raxml' else 2
        else:
            raise NotImplementedError
        self.save()
        sleep(.1)

        iface.thread = threading.Thread(target=self.do_ML, args=[mode])
        GObject.timeout_add(100, self.update_ML, PAGE, ml)
        iface.thread.start()
        return  # to main loop

    def do_ML(self, mode):
        """Run the ML inference thread"""
        data = self.data
        iface = self.iface
        ml = data.ml
        start = time()
        msa = self.wd / repo.PATHS.msa
        prefix = self.wd / 'RAxML'
        shell = self.wd / 'raxml_run.sh'
        errors = list()
        iface.i = 0

        # prepare the calls
        chck = '%s --msa %s --check --model %s' \
               + ml.evo_modify + ' --prefix %s'

        inML = '%s --msa %s --model %s' + ml.evo_modify + \
               ' --prefix %s' + ' --seed %d' % ml.raxml_seed + \
               ' --threads auto{16} --workers auto{16}' \
               ' --redo --tree %s ' % ','.join(
            [a for a in ['rand{%d}' % ml.rand if ml.rand > 0 else None,
                         'pars{%d}' % ml.pars if ml.pars > 0 else None] if a])

        boot = '%s --bootstrap --msa %s --model %s --tree %s' + \
               ' --prefix %s' + ' --bs-trees %d' % ml.bootstraps + \
               ' --seed %d' % ml.raxml_seed + \
               ' --threads auto{16} --workers auto{16} --redo'

        supp = '%s --support --tree %s --bs-trees %s --bs-metric fbp,tbe ' \
               '--prefix %s --threads auto{16} --workers auto{16} --redo'

        # res = run(stdout=PIPE, stderr=PIPE, args=shlex.split(
        #     arg % (ml.raxml, msa, prefix / 'bs_')))
        # res = run(stdout=PIPE, stderr=PIPE, shell=True,
        #           args=arg % (ml.raxml, msa, prefix / 'bs_'))
        # notify = 'notify-send "AB12PHYLO" "ML Tree Inference finished!" -u normal -i "%s"' \
        #          % str(BASE_DIR / 'ab12phylo' / 'files' / 'favi.png')
        # notify2 = 'zenity --notification --text="AB12PHYLO\nML Tree Inference finished" ' \
        #           '--window-icon="%s"' % str(BASE_DIR / 'ab12phylo' / 'files' / 'favi.png')

        if ml.raxml_shell:
            with open(shell, 'w') as sh:
                sh.write('#!/bin/bash\n\ntrap \'\' SIGINT\n\n')

        # loop over the stages
        for i, (desc, key, prev, arg, add) in enumerate(zip(
                ['check MSA', 'infer ML tree', 'bootstrapping', 'calc. branch support'],
                [False, 'ML', 'BS', False], [0, 1, ml.rand + ml.pars + 1, iface.k - 2],
                [chck, inML, boot, supp],
                [(ml.raxml, msa, ml.evo_model, prefix / 'chk'),
                 (ml.raxml, msa, ml.evo_model, prefix / 'ml'),
                 (ml.raxml, msa, prefix / 'ml.raxml.bestModel',
                  prefix / 'ml.raxml.bestTree', prefix / 'bs'),
                 (ml.raxml, prefix / 'ml.raxml.bestTree',
                  prefix / 'bs.raxml.bootstraps', prefix / 'sp')])):

            if ml.raxml_shell and mode == 'raxml':
                with open(shell, 'a') as sh:
                    sh.write('# %s\n' % desc)
                    sh.write(arg % add)
#                     if i != 2:
#                         sh.write(arg % add)
#                     else:
#                         # special case bootstrapping:
#                         bash = '''
# mkfifo pipe || exit 1
# (%s) > pipe &
# pid=$!
# echo "AB12PHYLO: Bootstrapping PID is $pid"
# while read -r line; do
#     echo "$line"
#     if [[ "${line::12}" == "Elapsed time" ]]; then
#         echo "AB12PHYLO: Finished bootstrapping, terminating process $pid to ensure it exits."
#         kill -s SIGTERM $pid
#         break
#     fi
# done < pipe
# rm pipe
#                         '''.strip() % (arg % add)
#                         sh.write(bash)
#                     sh.write('\n\necho "AB12PHYLO: %s done"\n' % desc)
                    if i != 3:
                        sh.write('\nsleep 1s\n\n')
                    continue

            iface.text = desc
            LOG.info(iface.text)
            ml.stdout = list()
            ml.key = key
            ml.prev = prev

            # read realtime RAxML output line by line
            proc = Popen(args=shlex.split(arg % add), stdout=PIPE, stderr=PIPE)
            while True:
                line = proc.stdout.readline()
                if proc.poll() is not None:
                    sleep(.2)
                    break
                if line:
                    lane = line.decode().rstrip()
                    ml.stdout.append(lane)
                    LOG.debug(lane)
                    if lane.startswith('Elapsed time'):
                        break
                else:
                    sleep(.2)
                    break

            # bf = iface.ml_help.get_buffer()
            # bf.props.text = bf.props.text + '\n' + '\n'.join(ml.stdout)
            # bf.insert_markup(bf.get_end_iter(),
            #                  '<span foreground="#2374AF">'
            #                  '________________________________________________'
            #                  '________________________________\n</span>', -1)
            # mark = bf.create_mark(None, bf.get_end_iter(), True)
            # iface.ml_help.scroll_mark_onscreen(mark)
            # bf.add_mark(Gtk.TextMark.new(stage, True), bf.get_end_iter())
            # iface.ml_help.scroll_to_iter(bf.get_end_iter(), .1, False, 0, .9)

            # check for errors
            for line in ml.stdout:
                if line.startswith('ERROR'):
                    errors.append(line)
            if errors:
                GObject.idle_add(self.stop_ML, errors, start)
                return True

            if mode == 'raxml_export':
                iface.text = 'building zip'
                LOG.debug(iface.text)
                sh = 'raxml_run.sh'
                with open(sh, 'w') as sf:
                    sf.write('#!/bin/bash')
                    sf.write('\n\n# Check MSA\n')
                    sf.write(chck % ('./raxml-ng', 'msa.fasta', Path(ml.evo_model).name, 'chk'))
                    sf.write('\n\n# Find best ML tree\n')
                    sf.write(inML % ('./raxml-ng', 'msa.fasta', Path(ml.evo_model).name, 'ml'))
                    sf.write('\n\n# Compute bootstrap iterations\n')
                    sf.write(boot % ('./raxml-ng', 'msa.fasta', 'ml.raxml.bestModel',
                                     'ml.raxml.bestTree', 'bs'))
                    sf.write('\n\n# Calculate branch support\n')
                    sf.write(supp % ('./raxml-ng', 'ml.raxml.bestTree',
                                     'bs.raxml.bootstraps', 'sp'))
                    sf.write('\n\n')

                with ZipFile(self.wd / 'RAxML_export.zip', 'w', ZIP_DEFLATED) as zf:
                    if ml.evo_model.endswith('user_model'):
                        zf.write(ml.evo_model, 'user_model')
                    zf.write(ml.raxml, 'raxml-ng')
                    zf.write(msa, 'msa.fasta')
                    zf.write(sh)
                Path(sh).unlink()

                sleep(.05)
                GObject.idle_add(self._save_somewhere_else, self.wd / 'RAxML_export.zip')
                return True

        if ml.raxml_shell and mode == 'raxml':
            # with open(shell, 'a') as sh:
            #     sh.write('# restart AB12PHYLO afterwards\n'
            #              'ab12phylo --open %s --proceed\n'
            #              % str(gui.wd / gui.project_path))
            shell.chmod(shell.stat().st_mode | stat.S_IEXEC)
            self.hold()
            self.win.hide()
            sleep(.1)
            Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference running in background.',
                   '-i', str(repo.PATHS.icon_path)])
            os.system(shell)
            self.win.show_all()
            self.release()
            # old alternative that killed python:
            # os.execv(shell, ('placeholder', 'arg'))

        iface.text = 'copy tree files'
        iface.i = iface.k - 1
        shutil.copy(prefix / 'sp.raxml.supportFBP', self.wd / repo.PATHS.fbp)
        shutil.copy(prefix / 'sp.raxml.supportTBE', self.wd / repo.PATHS.tbe)
        iface.text = 'idle'
        iface.frac = 1
        sleep(.1)
        GObject.idle_add(self.stop_ML, errors, start)
        return True

    def stop_ML(self, errors, start):
        """Finish the ML inference thread"""
        iface = self.iface
        iface.thread.join()
        self.update_ML(PAGE, self.data.ml)
        self.win.show_all()
        LOG.info('ML thread idle')
        self.set_changed(PAGE, False)
        if errors:
            self.show_notification('Errors during ML inference', errors)
        elif time() - start > 120:
            Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference finished',
                   '-i', str(repo.PATHS.icon_path)])
            # notify = threading.Thread(target=_zenity, args=())
            # notify.start()
        else:
            self.show_notification('ML finished')
        self.refresh()
        if iface.run_after:
            [do_func() for do_func in iface.run_after]
        return

    def _save_somewhere_else(self, path):
        dialog = Gtk.FileChooserDialog(title='export zip',
                                       parent=None, select_multiple=False,
                                       action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_filename(str(path))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            p = Path(dialog.get_filename()).resolve()
            try:
                shutil.move(path, p)
            except Exception as ex:
                LOG.error(ex)
        dialog.destroy()
