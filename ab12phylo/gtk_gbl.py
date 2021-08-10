# 2021 Leo Kaindl

import logging
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from time import sleep

import gi
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.cm import get_cmap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject

from ab12phylo import repo
from ab12phylo_cmd.filter import chmod_x
from ab12phylo.gtk_base import ab12phylo_app_base

LOG = logging.getLogger(__name__)
PAGE = 4


class gbl_page(ab12phylo_app_base):

    def __init__(self):
        super().__init__()
        data = self.data
        iface = self.iface

        iface.view_gbl.set_model(data.gbl_model)
        iface.view_gbl.append_column(Gtk.TreeViewColumn(
            title='id', cell_renderer=Gtk.CellRendererText(), text=0))
        iface.view_gbl.connect('check-resize', self.get_height_resize, iface.gbl_spacer,
                               [iface.gbl_left, iface.gbl_right])

        iface.gbl_preset.set_id_column(0)
        iface.gaps.set_id_column(0)

        iface.tempspace.bak_ignore = set()
        iface.tempspace.params = list()
        for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
            wi = iface.__getattribute__(w_name)
            iface.tempspace.__setattr__(w_name, wi.get_adjustment())
        for wi in [iface.conserved, iface.flank]:
            wi.connect_after('changed', self.edit)  # when editing field, check adjustments

        iface.gaps.connect('changed', lambda *args: data.gbl.  # no assignment in lambda -> use __setattr__
                           __setattr__('gaps', iface.gaps.get_active_text()))
        iface.gbl_preset.connect_after('changed', self.re_preset)
        iface.gbl_preset.set_active_id('balanced')

        sel = iface.view_gbl.get_selection()
        sel.set_mode(Gtk.SelectionMode.MULTIPLE)
        sel.connect('changed', self.keep_visible, iface.view_gbl,
                    iface.parallel_gbl.props.vadjustment.props, iface.tempspace)
        iface.view_gbl.connect('key_press_event', self.delete_and_ignore_rows,
                               PAGE, sel)  # in-preview deletion
        iface.gbl_eventbox.connect('button_press_event', self.select_seqs, PAGE, iface.zoomer,
                                   iface.view_gbl, iface.tempspace)  # in-preview selection
        iface.gbl_eventbox.connect('scroll-event', self.xy_scale, PAGE)  # zooming

    def refresh(self):
        """Cause the initial plotting or load .png images for good responsiveness later on"""
        data = self.data
        iface = self.iface

        if not data.gene_ids:
            return

        # start the first re-plotting
        if 0 in data.msa_shape or not (self.wd / repo.PATHS.left).exists() \
                or not (self.wd / repo.PATHS.right).exists():
            # fetch the MSA shape
            data.gbl_model.clear()
            shared_ids = sorted(set.intersection(*data.gene_ids.values()))
            [data.gbl_model.append([_id]) for _id in shared_ids]
            i = len(shared_ids)
            data.msa_shape = [sum(data.msa_lens), i, -1, i]  # width-height-width_after-height
            iface.msa_shape.set_text('%d : %d' % tuple(data.msa_shape[:2]))
            iface.msa_shape_trimmed.set_text('-1 : %d' % i)
            # set the adjustment boundaries of the spin buttons
            self.re_preset(iface.gbl_preset)
            self.start_gbl()
            return

        # else load the pre-existing PNGs
        try:
            iface.msa_shape.set_text('%d : %d' % tuple(data.msa_shape[:2]))
            iface.msa_shape_trimmed.set_text(
                '%d : %d' % (data.msa_shape[2] - (len(data.genes) - 1) * len(repo.SEP),
                             data.msa_shape[3]))
        except TypeError:
            pass
        # place the existing png
        x_ratio = data.msa_shape[2] / data.msa_shape[0]
        try:
            self.load_image(iface.zoomer, PAGE, iface.gbl_left_vp, self.wd / repo.PATHS.left,
                            w=data.gbl_shape[0] * self.get_hadj(),
                            h=data.gbl_shape[1])
            self.load_image(iface.zoomer, PAGE, iface.gbl_right_vp, self.wd / repo.PATHS.right,
                            w=data.gbl_shape[0] * self.get_hadj() * x_ratio,
                            h=data.gbl_shape[1])
            self.load_colorbar(iface.palplot2)
            self.re_preset(iface.gbl_preset)
            self.win.show_all()
        except ValueError:
            data.msa_shape[3] = 0
            self.refresh()

    def reload_ui_state(self):
        data = self.data
        iface = self.iface
        iface.gbl_preset.set_active_id(data.gbl.preset)
        iface.gaps.set_active_id(data.gbl.gaps)
        for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
            iface.tempspace.__getattribute__(w_name).configure(*data.gbl.__getattribute__(w_name), 1, 0, 0)
            # configure(value, lower, upper, step-increment=1, page-increment=0, page-size=0)

    def re_preset(self, gbl_preset):
        """
        Apply a preset on toggling it.
        :param gbl_preset: the GtkComboBoxText with the presets
        """
        mode = gbl_preset.get_active_text()
        g = self.iface.tempspace
        n_sites, n_seqs = self.data.msa_shape[:2]

        assert type(n_sites) == int
        # block / un-block SpinButtons
        [self.iface.__getattribute__(w_name).set_sensitive(mode != 'skip')
         for w_name in ['conserved', 'flank', 'good_block', 'bad_block', 'gaps']]

        if mode == 'skip':
            return
        elif mode == 'relaxed':
            conserved = n_seqs // 2 + 1
            flank = conserved
            gaps = 'half'
            good_block = 5
            bad_block = 8
        elif mode == 'balanced':
            conserved = n_seqs // 2 + 1
            flank = min(n_seqs // 4 * 3 + 1, n_seqs)
            gaps = 'half'
            good_block = 5
            bad_block = 8
        elif mode == 'default':
            conserved = n_seqs // 2 + 1
            flank = min(int(n_seqs * 0.85) + 1, n_seqs)
            gaps = 'none'
            good_block = 10
            bad_block = 8
        elif mode == 'strict':
            conserved = int(n_seqs * .9)
            flank = conserved
            gaps = 'none'
            good_block = 10
            bad_block = 8
        else:
            raise RuntimeWarning(f'illegal mode: {mode}')

        # configure(value, lower, upper, step-increment=1, page-increment=0, page-size=0)
        g.conserved.configure(conserved, n_seqs // 2 + 1, n_seqs, 1, 0, 0)
        g.flank.configure(flank, n_seqs // 2 + 1, n_seqs, 1, 0, 0)
        g.good_block.configure(good_block, 2, n_sites, 1, 0, 0)
        g.bad_block.configure(bad_block, 0, n_sites, 1, 0, 0)
        self.iface.gaps.set_active_id(gaps)

    def edit(self, wi):
        """Just make sure flank is never smaller than conserved"""
        LOG.debug('editing')
        adj = wi.get_adjustment()
        flank = self.iface.tempspace.flank
        conserved = self.iface.tempspace.conserved
        if adj == flank:
            return  # break some loops
        # conserved has changed
        flank.configure(max(adj.get_value(), flank.get_value()),  # value
                        adj.get_value(), flank.get_upper(), 1, 0, 0)  # the others

    def start_gbl(self, run_after=None, force=False):
        """Set-up the Gblocks thread. This cannot be reached if Gblocks shall be skipped"""
        data = self.data
        iface = self.iface
        if not data.genes:
            LOG.debug('abort Gblocks')
            return
        elif iface.thread.is_alive():
            self.show_notification('Busy', secs=1)
            return
        elif run_after and (self.wd / repo.PATHS.msa).exists() \
                and not self.get_errors(PAGE):
            self.set_changed(PAGE, False)
            [do_func() for do_func in run_after]
            return

        iface.refresh.grab_focus()  # trick to make the adjustments realize they've changed
        if iface.gbl_preset.get_active_text() != 'skip':
            self.set_changed(PAGE, True)
            # get arguments from adjustments values
            g = iface.tempspace  # the param Namespace
            # do NOT swap order, do NOT swap order! Gblocks wants them in the wrong order!
            flank, cons, good, bad = [adj.get_value() for adj in
                                      [g.flank, g.conserved, g.good_block, g.bad_block]]
            gaps = iface.gaps.get_active_text()[0]
            # failsafe
            if flank < cons:
                flank = cons
            params = [flank, cons, bad, good, gaps]
            # check if they have changed
            if iface.tempspace.params == params and not force \
                    and iface.tempspace.bak_ignore == data.gbl.ignore_ids \
                    and not self.get_changed(PAGE):
                self.show_notification('MSA already trimmed with these parameters', secs=2)
                return
            iface.tempspace.params = params
        data.gbl_model.clear()

        # save_ui_state from iface.tempspace into data.gbl
        data.gbl.preset = iface.gbl_preset.get_active_text()
        data.gbl.gaps = iface.gaps.get_active_text()
        for w_name in ['conserved', 'flank', 'good_block', 'bad_block']:
            p = iface.tempspace.__getattribute__(w_name).props
            data.gbl.__setattr__(w_name, [p.value, p.lower, p.upper])

        iface.thread = threading.Thread(target=self._do_gbl)
        iface.run_after = run_after
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()
        sleep(.1)
        # return to main loop

    def _do_gbl(self):
        """Run the Gblocks thread"""
        data = self.data
        iface = self.iface

        iface.frac = .05
        iface.i = 0
        iface.k = len(data.genes) + 6

        if iface.gbl_preset.get_active_text() != 'skip':
            iface.text = 'prep'
            LOG.debug(iface.text)
            # look for local system Gblocks
            binary = shutil.which('Gblocks')
            local = bool(binary)
            # else pick deployed Gblocks
            binary = Path(binary) if binary else next(repo.TOOLS.rglob(f'Gblocks{repo.EXE}'))
            # Make the file executable
            chmod_x(binary)
            LOG.info('%s Gblocks' % ('local' if local else 'packaged'))

            # create base call -t=d sets the mode to nucleotides ... adapt?
            # leave -b2 and -b1 in this wrong order!
            arg = '"%s" %s -t=d -b2=%d -b1=%d -b3=%d -b4=%d -b5=%s -e=".txt" -s=y -p=n -k=y; exit 0' \
                  % tuple([binary, '"%s"'] + iface.tempspace.params)
            if sys.platform in ['win32', 'cygwin']:
                arg += ' & exit /b 0'
            else:
                arg += '; exit 0'
            LOG.debug(arg)
        else:
            arg = ''

        # fetch IDs
        shared_ids = sorted(set.intersection(*data.gene_ids.values()) - data.gbl.ignore_ids)
        [data.gbl_model.append([_id]) for _id in shared_ids]

        iface.i += 1
        GObject.idle_add(self._do_gbl2, shared_ids, arg)
        sleep(.1)
        return True

    def _do_gbl2(self, shared_ids, arg):
        """Bounce to idle to size the id column"""
        data = self.data
        iface = self.iface
        iface.thread.join()

        gbar, barspace = (True, 50) if len(data.genes) > 1 else (False, 10)

        iface.view_gbl.set_margin_bottom(barspace)
        iface.view_gbl.realize()
        sleep(.1)
        data.gbl_shape[1] = self.get_height_resize(iface.view_gbl, None, iface.gbl_spacer,
                                                   [iface.gbl_left, iface.gbl_right]) + barspace

        # re-direct to thread
        iface.thread = threading.Thread(target=self._do_gbl3, args=[shared_ids, arg, gbar, barspace])
        GObject.timeout_add(100, self.update, PAGE)
        iface.thread.start()

    def _do_gbl3(self, shared_ids, arg, gbar, barspace):
        data = self.data
        iface = self.iface

        errors = list()
        msa_lens = list()
        blocks = list()
        array = np.empty(shape=(len(shared_ids), 0), dtype=int)

        for gene in data.genes:
            iface.text = '%s: read MSA' % gene
            LOG.debug(iface.text)
            raw_msa = (self.wd / gene / ('%s_raw_msa.fasta' % gene)).resolve()
            msa = self.wd / gene / ('%s_msa.fasta' % gene)
            records = {r.id: r for r in SeqIO.parse(raw_msa, 'fasta')}
            take_out = {_id for _id in records.keys() if _id in data.gbl.ignore_ids}
            take_out = {_id: records.pop(_id) for _id in take_out}
            if take_out:
                # write newly dropped sequences to backup file
                new_take_out = {_id: r for _id, r in take_out.items()
                                if _id not in iface.tempspace.bak_ignore}
                with open(self.wd / gene / ('%s_raw_msa_dropped.fasta' % gene), 'a') as fasta:
                    SeqIO.write(new_take_out.values(), fasta, 'fasta')
                # overwrite MSA without all dropped samples
                with open(raw_msa, 'w') as fasta:
                    SeqIO.write(records.values(), fasta, 'fasta')
            if sorted(records.keys()) != shared_ids:
                errors.append('MSA for %s does not match the dataset, please re-build.' % gene)
                sleep(.1)
                GObject.idle_add(self.stop_gbl, errors)
                return True
            ar = np.array([repo.seqtoint(records[_id].seq.upper()) for _id in shared_ids])
            msa_lens.append(ar.shape[1])
            # get the array columns that are not only gaps
            usable_sites = [i for i in range(ar.shape[1]) if set(ar[:, i]) != {repo.toint('-')}]
            array = np.hstack((array, ar,))  # shared.SEP would need to be stacked here
            del ar
            data.msa_shape[:2] = array.shape[::-1]

            if iface.gbl_preset.get_active_text() == 'skip':
                LOG.debug('skipping %s' % gene)
                shutil.copy(raw_msa, msa)
                iface.i += 1
                continue

            iface.text = '%s: run Gblocks' % gene
            LOG.debug(iface.text)
            with open(self.wd / gene / 'gblocks.log', 'w') as log_handle:
                try:
                    LOG.debug(arg % raw_msa)
                    subprocess.run(arg % raw_msa, shell=True, check=True,
                                   stdout=log_handle, stderr=log_handle)
                except (OSError, subprocess.CalledProcessError) as e:
                    errors.append(str(e))
                    log_handle.write(str(e))
                    continue

            # parse result
            iface.text = '%s: parse result' % gene
            LOG.debug(iface.text)

            shutil.move(raw_msa.with_suffix('.fasta.txt'), msa)
            # get the good blocks from the last pseudo-sequence in the text mask file
            for pseudo_seq in SeqIO.parse(raw_msa.with_suffix('.fasta.txtMask'), 'fasta'):
                mask = pseudo_seq
            # map the Gblocks mask to the original MSA sites
            line_blocks = [usable_sites[i] for i, char in enumerate(mask.seq) if char == '#']
            # LOG.debug(line_blocks)
            if not line_blocks:
                err = '%s: no good blocks' % gene
                LOG.error(err)
                errors.append(err)
                continue

            shift = sum(msa_lens[:-1])
            blocks.extend([i + shift for i in line_blocks])
            iface.i += 1
        data.msa_lens = msa_lens

        iface.text = 'concatenating MSAs'
        LOG.debug(iface.text)
        iface.tempspace.bak_ignore = {i for i in data.gbl.ignore_ids}
        # if 'aligner' not in iface:
        iface.aligner, cmd = self.get_msa_build_cmd(
            repo.toalgo(iface.msa_algo.get_active_text()), self.wd, data.genes)
        # iface.aligner.reset_paths(self.wd, self.wd / repo.PATHS.msa)
        data.msa_shape[2], data.msa_shape[3] = iface.aligner.concat_msa(gui=shared_ids)
        iface.text = 'computing SHA256 hash'
        LOG.debug(iface.text)
        self.get_hashes(repo.PATHS.msa, PAGE)

        iface.i += 1
        iface.text = 'plot MSAs'
        LOG.debug(iface.text)

        data.gbl_shape[0] = array.shape[1] * self.get_hadj()
        # make gaps transparent
        array = np.ma.masked_where(array > repo.toint('else'), array)

        # create a transparency mask
        gbl_mask = np.full(array.shape, repo.ALPHA)
        gbl_mask[:, blocks] = 1

        data.msa_shape[2] = len(blocks)  # for completeness
        LOG.debug('msa shape: %s' % str(data.msa_shape))
        x_ratio = data.msa_shape[2] / data.msa_shape[0]
        LOG.debug('x ratio: %.3f' % x_ratio)

        # adjust maximum size
        scale = 6
        while max(data.msa_shape[:2]) * scale > 2 ** 14:
            scale -= 1
        LOG.debug('scaling gbl with %d' % scale)

        if iface.gbl_preset.get_active_text() != 'skip':
            for alpha, blocks, gtk_bin, png_path, x_ratio, width \
                    in zip([gbl_mask, 1], [range(array.shape[1]), blocks],
                           [iface.gbl_left_vp, iface.gbl_right_vp],
                           [repo.PATHS.left, repo.PATHS.right], [1, x_ratio],
                           [data.msa_shape[0], data.msa_shape[2]]):
                f = Figure()  # figsize=(width / shared.DPI, data.msa_shape[1] / shared.DPI), dpi=shared.DPI)  # figaspect(data.msa_shape[1] / width))
                f.set_facecolor('none')
                f.set_figheight(data.msa_shape[1] / repo.DPI * 5)
                f.set_figwidth(max(1, width) / repo.DPI)
                # f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

                # leave room at the bottom for the gene legend bar
                b = barspace / data.gbl_shape[1]
                ax = f.add_axes([0, b, 1, 1 - b])

                mat = ax.matshow(array[:, blocks], alpha=alpha, cmap=ListedColormap(repo.colors),
                                 vmin=-.5, vmax=len(repo.colors) - .5, aspect='auto')
                if not gbar:
                    LOG.debug('adding ticks')
                    [ax.spines[t].set_visible(False) for t in ['left', 'right', 'top', 'bottom']]
                    ax.yaxis.set_visible(False)
                    ax.xaxis.set_ticks_position('bottom')
                    ax.tick_params(colors=iface.FG, pad=.2, length=0, labelsize=1)
                    ax.xaxis.set_ticks([i for i in range(1, array[:, blocks].shape[1] - 1) if not i % 100])
                else:
                    ax.axis('off')
                    LOG.debug('adding gene marker bar')
                    # build matrix of gene indicators
                    gm = [[i] * l for i, l in enumerate(msa_lens)]
                    if x_ratio == 1:
                        # add the spacer
                        for i in range(len(data.genes) - 1):
                            gm[i] += [-1] * len(repo.SEP)
                    gm = [gi for gl in gm for gi in gl]
                    # turn into np array
                    gm = np.vstack([np.array(gm)] * 2)
                    # make spacer transparent
                    gm = np.ma.masked_where(gm < 0, gm)
                    # trim the array early
                    gm = gm[:, blocks]
                    gene_colors = get_cmap('GnBu', len(data.genes))

                    # plot the marker bar onto the MSA graphic
                    bax = f.add_axes([0, b * .4, 1, b / 3])
                    bar = bax.pcolormesh(gm, cmap=gene_colors)
                    # bax.axis('off')
                    [bax.spines[t].set_visible(False) for t in ['left', 'right', 'top', 'bottom']]
                    bax.yaxis.set_visible(False)
                    bax.xaxis.set_ticks_position('bottom')
                    bax.tick_params(colors=iface.FG, pad=.2, length=0, labelsize=1)
                    bax.xaxis.set_ticks([i for i in range(1, gm.shape[1] - 1) if not i % 100])

                    # plot a legend for the gene marker bar
                    with plt.rc_context({'axes.edgecolor': iface.FG, 'xtick.color': iface.FG}):
                        iface.text = 'gene marker bar'
                        LOG.debug(iface.text)
                        Path.mkdir(self.wd / repo.PATHS.phylo_msa.parent, exist_ok=True)
                        fig = plt.figure(figsize=(4, .2))
                        cax = fig.add_subplot(111)
                        cbar = ColorbarBase(
                            ax=cax, cmap=gene_colors, orientation='horizontal',
                            ticks=[(.5 / len(data.genes) + j * 1 / len(data.genes))
                                   for j in range(len(data.genes))])
                        cbar.ax.set_xticklabels(data.genes)
                        fig.savefig(self.wd / repo.PATHS.gbar, transparent=True,
                                    bbox_inches='tight', pad_inches=0, dpi=600)
                        plt.close(fig)
                        del fig, cbar

                iface.i += 1
                iface.text = 'save PNG'
                LOG.debug(iface.text)
                Path.mkdir(self.wd / png_path.parent, exist_ok=True)
                f.savefig(self.wd / png_path, transparent=True,
                          dpi=scale * repo.DPI, bbox_inches='tight', pad_inches=0.00001)

                if iface.rasterize.props.active:
                    iface.text = 'place PNG'
                    LOG.debug(iface.text)
                    self.load_image(iface.zoomer, PAGE, gtk_bin, self.wd / png_path,
                                    data.gbl_shape[0] * x_ratio * self.get_hadj(), data.gbl_shape[1])
                else:
                    iface.text = 'place vector'
                    LOG.debug(iface.text)
                    canvas = FigureCanvas(f)
                    canvas.set_size_request(max(len(blocks) * self.get_hadj(), -1),
                                            data.gbl_shape[1])  # width, height
                    try:
                        ch = gtk_bin.get_child()
                        if ch:
                            gtk_bin.remove(ch)
                        gtk_bin.add(canvas)
                    except Exception as ex:
                        LOG.error(ex)
                iface.i += 1

                gtk_bin.realize()
                gtk_bin.show_all()

        # re-size
        LOG.debug('re-sizing again')
        for wi in [iface.gbl_left, iface.gbl_right]:
            wi.set_max_content_height(data.gbl_shape[1])

        self.load_colorbar(iface.palplot2)

        if gbar:
            self.load_colorbar(iface.gbar1, gbar=True)
            iface.gbar1.set_visible(True)
        else:
            iface.gbar1.set_visible(False)

        iface.text = 'idle'
        iface.frac = 1
        sleep(.1)
        GObject.idle_add(self.stop_gbl, errors)
        return True

    def stop_gbl(self, errors):
        """Finish the Gblocks thread"""
        data, iface = self.data, self.iface
        iface.thread.join()
        LOG.info('gbl thread idle')
        sleep(.1)
        self.update(PAGE)
        self.win.show_all()
        self.set_errors(PAGE, bool(errors))
        self.set_changed(PAGE, False)
        iface.msa_shape.set_text('%d : %d' % tuple(data.msa_shape[:2]))
        iface.msa_shape_trimmed.set_text(
            '%d : %d' % (data.msa_shape[2], data.msa_shape[3]))

        # possibly re-configure flank and conserved
        g = iface.tempspace
        n_seqs = data.msa_shape[1]
        min_c = n_seqs // 2 + 1
        if g.flank.props.upper != n_seqs:
            g.conserved.configure(max(min_c, g.conserved.props.value), min_c, n_seqs, 1, 0, 0)
            g.flank.configure(max(min_c, g.flank.props.value), min_c, n_seqs, 1, 0, 0)
        if errors:
            self.show_notification('Errors during MSA trimming', errors)
            return
        if iface.run_after:
            [do_func() for do_func in iface.run_after]
        return

    def undrop_seqs(self):
        """Write dropped sequences back to the untrimmed MSAs."""
        ignore_ids, data, iface = self.data.gbl.ignore_ids, self.data, self.iface
        for gene in data.genes:
            dropped_file = self.wd / gene / ('%s_raw_msa_dropped.fasta' % gene)
            if not dropped_file.is_file():
                continue
            put_in = {r.id: r for r in SeqIO.parse(dropped_file, 'fasta')}

            unknown = {_id for _id in put_in.keys() if _id not in ignore_ids}
            assert not unknown, 'Found wrong dropped sequences: %s' % ','.join(unknown)

            with open(self.wd / gene / ('%s_raw_msa.fasta' % gene), 'a') as fasta:
                SeqIO.write(put_in.values(), fasta, 'fasta')
            dropped_file.unlink()

        # save in metadata
        for _id in ignore_ids:
            for gene, genemeta in data.metadata.items():
                if _id in genemeta:
                    genemeta[_id]['quality'] = ''
        self.write_metadata()

        ignore_ids.clear()
        self.iface.tempspace.bak_ignore.clear()
        self.start_gbl(force=True)

    def drop_seqs(self):
        """
        On manually deleting samples from the trimmed MSA, overwrite 'old'
        trimmed MSA file, save in and write metadata.
        """
        ignore_ids, data, iface = self.data.gbl.ignore_ids, self.data, self.iface
        if not ignore_ids:
            return
        # overwrite old MSA
        with open(self.wd / (repo.PATHS.msa + '_TEMP'), 'w') as fasta:
            for record in SeqIO.parse(self.wd / repo.PATHS.msa, 'fasta'):
                if record.id not in ignore_ids:
                    SeqIO.write(record, fasta, 'fasta')
        LOG.debug('dropped %d sequences from trimmed MSA' % len(ignore_ids))
        shutil.move(self.wd / (repo.PATHS.msa + '_TEMP'), self.wd / repo.PATHS.msa)

        # save in metadata
        for _id in ignore_ids:
            for gene, genemeta in data.metadata.items():
                if _id in genemeta:
                    genemeta[_id]['quality'] = 'manually dropped at trim_MSA'
        self.write_metadata()
