# 2020 Leo Kaindl

import logging
import os
import random
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

from ab12phylo_gui import shared
from ab12phylo_gui.static import PATHS as pp, TOOLS

# from ab12phylo_gui.gtk_main import ab12phylo_app

BASE_DIR = Path(__file__).resolve().parents[2]
LOG = logging.getLogger(__name__)
PAGE = 6


def init(gui):
    """Initialize the page. Connect buttons"""
    data, iface = gui.data, gui.iface
    iface.evo_model.set_model(data.evo_models)
    iface.evo_model.get_child().connect(
        'focus_out_event', _change_evo_model, iface.evo_modify, data.ml)
    iface.evo_block = iface.evo_model.connect_after(
        'changed', _load_model_file, iface.evo_modify, gui)
    iface.evo_model.handler_block(iface.evo_block)

    iface.raxml_run.connect('clicked', start_ML, gui, 'raxml')
    iface.raxml_export.connect('clicked', start_ML, gui, 'raxml_export')
    iface.raxml_import.connect('clicked', import_tree, gui)

    for wi in [iface.bootstraps, iface.rand, iface.pars, iface.raxml_seed]:
        wi.connect('key-press-event', shared.edit_numerical_entry_up_down)
        wi.connect('changed', shared.edit_numerical_entry)

    iface.raxml_seen = False


def import_tree(widget, gui):
    """Import a tree file or two"""
    data, ml, iface = gui.data, gui.data.ml, gui.iface

    dialog = Gtk.FileChooserDialog(title='import tree(s)',
                                   parent=None, select_multiple=True,
                                   action=Gtk.FileChooserAction.OPEN)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                       Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
    response = dialog.run()
    if response == Gtk.ResponseType.OK:
        paths = [Path(p).resolve() for p in dialog.get_filenames()]
        if len(paths) > 2:
            shared.show_message_dialog('Please select at most two tree files.')
        else:
            errors = list()
            for path in paths:
                if 'FBP' in path.name.upper():
                    shutil.copy(path, gui.wd / p.fbp)
                elif 'TBE' in path.name.upper():
                    shutil.copy(path, gui.wd / p.tbe)
                else:
                    errors.append(path.name)
            if errors:
                shared.show_message_dialog('Not immediately recognized as either '
                                           'tree_FBP.nwk or tree_TBE.nwk. You can also copy '
                                           'it/them to %s manually.' % gui.wd, list_to_print=errors)
            else:
                shared.show_notification(gui, 'imported trees:',
                                         [path.name for path in paths], 2)

    dialog.destroy()
    shared.set_changed(gui, PAGE, False)


def _change_evo_model(entry, event_focus, evo_modify, ml):
    combo = entry.get_parent().get_parent()
    if combo.get_active_iter():
        ml.evo_model = entry.get_text()
        evo_modify.props.sensitive = True
    else:
        tx = entry.get_text()
        combo.get_model().append([tx, None])
        ml.evo_model = tx  # save in project dataset
        LOG.debug('entered custom evo model %s' % tx)


def _load_model_file(combo, evo_modify, gui):
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
                Path.mkdir(gui.wd / 'RAxML', exist_ok=True)
                shutil.copy(tx, gui.wd / 'RAxML' / 'user_model')
            except FileNotFoundError:
                shared.show_notification(gui, 'File does not exist', stay_secs=1)
                dialog.destroy()
                return
            except shutil.SameFileError:
                pass
            evo_modify.props.sensitive = False
            combo.get_model().append([tx.name, str(tx)])
            combo.set_active_id(tx.name)
            LOG.debug('selected partitioned model file: %s -> RAxML/user_model' % str(tx))
        dialog.destroy()


def _load_raxml_help(gui):
    data, iface = gui.data, gui.iface
    binary = shutil.which('raxml-ng')
    if not binary:
        binary = str(TOOLS / pp.RAxML)
    iface.raxml_exe.set_filename(binary)
    res = run(stdout=PIPE, stderr=PIPE, shell=True,
              args='%s --help' % binary)
    iface.ml_help.props.buffer.props.text = res.stdout.decode().lstrip()
    LOG.debug('got RAxML --help')


def reload_ui_state(gui):
    data, ml, iface = gui.data, gui.data.ml, gui.iface
    for w_name in ['bootstraps', 'rand', 'pars']:
        iface.__getattribute__(w_name).set_text(str(ml.__getattribute__(w_name)))
    iface.raxml_seed.set_text(str(ml.raxml_seed) if 'raxml_seed' in ml else '')
    iface.evo_model.set_active_id(ml.evo_model)
    iface.evo_modify.set_text(ml.evo_modify)
    iface.raxml_shell.set_active(ml.raxml_shell)


def refresh(gui):
    """Re-view the page. Get suggested commands for RAxML-NG and IQ-Tree"""
    LOG.debug('ML refresh')
    data, iface = gui.data, gui.iface
    if not iface.raxml_seen:
        _load_raxml_help(gui)
        reload_ui_state(gui)
        iface.evo_model.handler_unblock(iface.evo_block)
        iface.evo_model.set_active_id(data.ml.evo_model)
        iface.raxml_seen = True

    shared.set_errors(gui, PAGE, not any((gui.wd / a).is_file() for a in [pp.tbe, pp.fbp]))


def start_ML(widget, gui, mode, run_after=None):
    """Set-up the ML inference thread"""
    data, ml, iface = gui.data, gui.data.ml, gui.iface

    if not data.genes or iface.thread.is_alive():
        shared.show_notification(gui, 'Busy', stay_secs=1)
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
            ml.evo_model = str(gui.wd / 'RAxML' / 'user_model')
        else:
            ml.evo_model = ml.evo_model[0]

        ml.raxml_seed = random.randint(0, max(1000, ml.bootstraps)) \
            if ml.raxml_seed == '' else int(ml.raxml_seed)
        iface.raxml_seed.props.text = str(ml.raxml_seed)
        Path.mkdir(gui.wd / 'RAxML', exist_ok=True)
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
    # ab12phylo_app.save(gui)  # TODO solve after refactoring

    iface.thread = threading.Thread(target=do_ML, args=[gui, mode])
    GObject.timeout_add(100, shared.update_ML, iface, PAGE, ml)
    iface.thread.start()
    return  # to main loop


def do_ML(gui, mode):
    """Run the ML inference thread"""
    start = time()
    data, ml, iface = gui.data, gui.data.ml, gui.iface
    msa = gui.wd / pp.msa
    prefix = gui.wd / 'RAxML'
    shell = gui.wd / 'raxml_run.sh'
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
    #          % str(BASE_DIR / 'ab12phylo_gui' / 'files' / 'favi.png')
    # notify2 = 'zenity --notification --text="AB12PHYLO\nML Tree Inference finished" ' \
    #           '--window-icon="%s"' % str(BASE_DIR / 'ab12phylo_gui' / 'files' / 'favi.png')

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
                if i != 2:
                    sh.write(arg % add)
                else:
                    # special case bootstrapping:
                    bash = '''     
mkfifo pipe || exit 1
(%s) > pipe &
pid=$!
echo "AB12PHYLO: Bootstrapping PID is $pid"
while read -r line; do
    echo "$line"
    if [[ "${line::12}" == "Elapsed time" ]]; then
        echo "AB12PHYLO: Finished bootstrapping, terminating process $pid to ensure it exits."
        kill -s SIGTERM $pid
        break
    fi
done < pipe
rm pipe
                    '''.strip() % (arg % add)
                    sh.write(bash)
                sh.write('\n\necho "AB12PHYLO: %s done"\n' % desc)
                if i != 3:
                    sh.write('sleep 1s\n\n')
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
            GObject.idle_add(stop_ML, gui, errors, start)
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

            with ZipFile(gui.wd / 'RAxML_export.zip', 'w', ZIP_DEFLATED) as zf:
                if ml.evo_model.endswith('user_model'):
                    zf.write(ml.evo_model, 'user_model')
                zf.write(ml.raxml, 'raxml-ng')
                zf.write(msa, 'msa.fasta')
                zf.write(sh)
            Path(sh).unlink()

            sleep(.05)
            GObject.idle_add(_save_somewhere_else, gui.wd / 'RAxML_export.zip')
            return True

    if ml.raxml_shell and mode == 'raxml':
        # with open(shell, 'a') as sh:
        #     sh.write('# restart AB12PHYLO afterwards\n'
        #              'ab12phylo-gui --open %s --proceed\n'
        #              % str(gui.wd / gui.project_path))
        shell.chmod(shell.stat().st_mode | stat.S_IEXEC)
        gui.hold()
        gui.win.hide()
        sleep(.1)
        Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference running in background.',
               '-i', str(Path(__file__).resolve().parent / 'files' / 'favi.png')])
        os.system(shell)
        gui.win.show_all()
        gui.release()
        # old alternative that killed python:
        # os.execv(shell, ('placeholder', 'arg'))

    iface.text = 'copy tree files'
    iface.i = iface.k - 1
    shutil.copy(prefix / 'sp.raxml.supportFBP', gui.wd / pp.fbp)
    shutil.copy(prefix / 'sp.raxml.supportTBE', gui.wd / pp.tbe)
    iface.text = 'idle'
    iface.frac = 1
    sleep(.1)
    GObject.idle_add(stop_ML, gui, errors, start)
    return True


def stop_ML(gui, errors, start):
    """Finish the ML inference thread"""
    iface = gui.iface
    iface.thread.join()
    shared.update_ML(iface, PAGE, gui.data.ml)
    gui.win.show_all()
    LOG.info('ML thread idle')
    shared.set_changed(gui, PAGE, False)
    if errors:
        shared.show_notification(gui, 'Errors during ML inference', errors)
    elif time() - start > 120:
        Popen(['notify-send', 'AB12PHYLO', 'ML Tree Inference finished',
               '-i', str(Path(__file__).resolve().parent / 'files' / 'favi.png')])
        # notify = threading.Thread(target=_zenity, args=())
        # notify.start()
    else:
        shared.show_notification(gui, 'ML finished')
    refresh(gui)
    if iface.run_after:
        [do_func(gui) for do_func in iface.run_after]
    return


def _save_somewhere_else(path):
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
