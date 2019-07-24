#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 11:10:47 2019

@author: Lau Møller Andersen
for chapter for Electrifying Language
tested used mne version '0.18.2
"""

#%% ===========================================================================
# SET PATHS AND IMPORTS
# =============================================================================

from os.path import join
import mne
import matplotlib.pyplot as plt
import matplotlib as mpl
from mayavi import mlab


path = '/home/lau/analyses/fieldtrip_datasets/language_MEG/'
data_path = join(path, 'data')
figures_path = join(path, 'figures')
subjects_dir = join(data_path, 'freesurfer')

subject = 'sub-V1002'
subject_path = join(data_path, subject, 'meg', subject + '_task-visual_meg.ds')
bem_folder = join(subjects_dir, subject, 'bem')
label_folder = join(subjects_dir, subject, 'label')
save_path = join(data_path, subject, 'meg', 'derived')

#%% ===========================================================================
# FIND EVENTS IN RAW DATA
# =============================================================================

raw = mne.io.read_raw_ctf(subject_path, preload=False)

events = mne.find_events(raw)
## there was a known delay between trigger and stimulation of 36 ms
shifted_events = mne.event.shift_time_events(events, [2, 4, 6, 8], 0.036,
                                             raw.info['sfreq'])
event_id = dict(rc_plus_sentence=2, rc_plus_scrambled=4,
                rc_minus_sentence=6, rc_minus_scrambled=8)

events = mne.write_events(join(save_path, 'reading-eve.fif'),
                          shifted_events)

#%% ===========================================================================
# EPOCH THE RAW DATA
# =============================================================================

tmin = -0.150
tmax =  0.500
baseline = (None, 0)
picks = mne.pick_types(raw.info, meg=True)

epochs = mne.Epochs(raw, shifted_events, event_id, tmin, tmax, baseline, picks)

mne.epochs.combine_event_ids(epochs,
                             ['rc_minus_sentence', 'rc_plus_sentence'],
                             dict(sentence=16), copy=False)
mne.epochs.combine_event_ids(epochs,
                             ['rc_minus_scrambled', 'rc_plus_scrambled'],
                             dict(scrambled=32), copy=False)

epochs.save(join(save_path, 'reading-epo.fif'))

#%% ===========================================================================
# EVOKEDS
# =============================================================================

evokeds = []
for event in epochs.event_id:
    evokeds.append(epochs[event].average())
    
mne.write_evokeds(join(save_path, 'reading-ave.fif'), evokeds)    
    
#%% ===========================================================================
# GRAPHICAL DEFAULTS
# =============================================================================

mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['font.size'] = 14
    
#%% ===========================================================================
# PLOT EVOKEDS
# =============================================================================

evokeds = mne.read_evokeds(join(save_path, 'reading-ave.fif'))

plt.close('all')
n_plots = 2
fig = plt.figure(figsize=(16, 9))
for plot_index in range(n_plots):
    plt.subplot(1, 2, plot_index + 1)
axes = fig.axes

for evoked_index, evoked in enumerate(evokeds):
    evoked.plot(ylim=dict(mag=(-300, 300)), axes=axes[evoked_index])
plt.savefig(join(figures_path, 'evoked.jpeg'), dpi=300)

mne.viz.plot_evoked_topomap

for evoked in evokeds:
    fig = evoked.plot_topomap(times=[0.100, 0.400], vmin=-120, vmax=120,
                              colorbar=True,
                              title=evoked.comment.capitalize())   
    fig.set_size_inches(9, 4)
    cbar = fig.axes[-1]
    text = cbar.set_ylabel('Magnetic Field Strength')
    tick_fonts = cbar.get_yticklabels()
    for tick_font in tick_fonts:
        tick_font.set_fontsize(14)

    plt.savefig(join(figures_path, 'evoked_topo_' + evoked.comment + '.jpeg'),
                dpi=300)

#%% ===========================================================================
# MR PREPROCESSING    
# =============================================================================

## watershed 
    
mne.bem.make_watershed_bem(subject, subjects_dir)

## source space

spacing = 'oct6'    
filename = subject + '-' + spacing[:3] + '-' + spacing[-1] + '-src.fif'
fullpath = join(bem_folder, filename)

src = mne.source_space.setup_source_space(subject, spacing,
                                          subjects_dir=subjects_dir,
                                          n_jobs=-1, surface='white')
mne.source_space.write_source_spaces(fullpath, src)

## bem surfaces
ico = 4
if ico == 3:
    ico_string = '1280'
elif ico == 4:
    ico_string = '5120'
elif ico == 5:
    ico_string = '20484'
filename = subject + '-' + ico_string + '-bem.fif'
fullpath = join(bem_folder, filename)

surfaces = mne.bem.make_bem_model(subject, ico, conductivity=[0.3],
                                      subjects_dir=subjects_dir)
mne.bem.write_bem_surfaces(fullpath, surfaces)
   
## bem solution
filename = subject + '-' + ico_string + '-bem-sol.fif'
fullpath = join(bem_folder, filename)
    
bem = mne.bem.make_bem_solution(surfaces)
mne.bem.write_bem_solution(fullpath, bem)

#%% ===========================================================================
# PLOT SOURCE SPACE
# =============================================================================

mlab.close(scene=None, all=True)

fig = mne.SourceSpaces.plot(src, subjects_dir=subjects_dir,
                            skull='inner_skull')
mlab.figure(fig, bgcolor=(1, 1, 1))
fig.scene.camera.zoom(1.5)
fig.scene.x_minus_view()
filename = 'bem_source_space_x_minus.jpeg'
fig.scene.save(join(figures_path, filename))
fig.scene.x_plus_view()
filename = 'bem_source_space_x_plus.jpeg'
fig.scene.save(join(figures_path, filename))
fig.scene.z_plus_view()
filename = 'bem_source_space_z_plus.jpeg'
fig.scene.save(join(figures_path, filename))
fig.scene.z_minus_view()
filename = 'bem_source_space_z_minus.jpeg'
fig.scene.save(join(figures_path, filename))


#%% ===========================================================================
# FORWARD SOLUTION
# =============================================================================
          
info = raw.info
trans = join(bem_folder, subject + '-trans.fif')
src = join(bem_folder, subject + '-' + spacing[:3] + '-' + spacing[-1] + \
           '-src.fif')
bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')
meg = True
eeg = False

filename = 'reading-fwd.fif'
fullpath = join(bem_folder, filename)
    
fwd = mne.make_forward_solution(info, trans, src, bem, meg, eeg,
                                n_jobs=-1)
mne.write_forward_solution(fullpath, fwd)

#%% ===========================================================================
# INVERSE OPERATOR
# =============================================================================

epochs = mne.read_epochs(join(save_path, 'reading-epo.fif'))

cov = mne.compute_covariance(epochs, tmax=0)
cov.save(join(save_path, 'reading-cov.fif'))

inv = mne.minimum_norm.make_inverse_operator(raw.info, fwd, cov)

mne.minimum_norm.write_inverse_operator(join(save_path, 'reading-inv.fif'),
                                        inv)

#%% ===========================================================================
# MNE
# =============================================================================

evokeds = mne.read_evokeds(join(save_path, 'reading-ave.fif'))
inv = mne.minimum_norm.read_inverse_operator(join(save_path,
                                                  'reading-inv.fif'))
stcs = dict()

for evoked in evokeds:
    stcs[evoked.comment] = mne.minimum_norm.apply_inverse(evoked, inv,
                                                              method='MNE')

#%% ===========================================================================
# PLOT MNE
# =============================================================================

mlab.close(None, True)

brain = mne.viz.plot_source_estimates(stcs['sentence'], subject, hemi='both',
                                      subjects_dir=subjects_dir,
                                      initial_time=0.100, background='w',
                                      foreground='k',
                                      clim=dict(kind='value',
                                                lims=[5e-11, 7e-11, 10e-11]),
                                      colorbar=True)
brain.show_view('caudal')
brain.save_image(join(figures_path, 'MNE_V1_sentence.jpeg'))

mlab.view(90, 135)
brain.set_time(0.400)
brain.save_image(join(figures_path, 'MNE_LOT_sentence.jpeg'))

#%% ===========================================================================
# LABELS
# =============================================================================

hemis = ['lh', 'rh']
labels = dict()
ltcs = dict()

for hemi in hemis:
    labels[hemi] = dict()
    labels[hemi]['V1'] = mne.read_label(join(label_folder,
          hemi + '.V1.label'))
    labels[hemi]['LOF'] = mne.read_label(join(label_folder,
          hemi + '.lateralorbitofrontal.label'))
    
ltcs = dict()
for hemi in hemis:
    ltcs[hemi] = dict()
    for label_index, label in enumerate(labels[hemi]):
        ltcs[hemi][label] = dict()
        for event_index, event in enumerate(stcs):
            ltcs[hemi][label][event] = \
                                stcs[event].in_label(labels[hemi][label])
                                
                                
#%% ===========================================================================
# PLOT LTCS                                
# =============================================================================
 
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['figure.titlesize'] = 30
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
    
fig = plt.figure(figsize=(16, 12))
for hemi_index, hemi in enumerate(ltcs):
    for label_index, label in enumerate(ltcs[hemi]):
        plot_number = 2*label_index + hemi_index + 1
        plt.subplot(2, 2, plot_number)
        times = stcs['sentence'].times * 1e3
        data = ltcs[hemi][label]['sentence'].data
        data = mne.filter.filter_data(data, 1200.0, None, 45) * 1e9
        plt.plot(times, data.mean(axis=0))
        data = ltcs[hemi][label]['scrambled'].data
        data = mne.filter.filter_data(data, 1200.0, None, 45) * 1e9
        plt.plot(times, data.mean(axis=0))
        plt.xlim()
        plt.ylim((0, 0.1))
        if plot_number > 2:
            plt.xlabel('Time (ms)')
        plt.ylabel('Current (Am)')
        plt.vlines(0, 0, 0.1, linestyles='dashed')
        plt.title(hemi.upper() + ' ' + label)
        if plot_number == 1:
            plt.legend(['Sentence', 'Scrambled'])
plt.suptitle('Average activity over ROIs')
plt.savefig(join(figures_path, 'label_time_courses.jpeg'), dpi=300)
                                
#%% ===========================================================================
# MNE for dSPM - not covered in tutorial
# =============================================================================

evokeds = mne.read_evokeds(join(save_path, 'reading-ave.fif'))
inv = mne.minimum_norm.read_inverse_operator(join(save_path,
                                                  'reading-inv.fif'))
stcs_dSPM = dict()
stcs_MNE = dict()

for evoked in evokeds:
    stcs_dSPM[evoked.comment] = mne.minimum_norm.apply_inverse(evoked, inv,
                                                              method='dSPM')
    stcs_MNE[evoked.comment] = mne.minimum_norm.apply_inverse(evoked, inv,
                                                              method='MNE')
    
stcs_diff_dSPM = stcs_dSPM['sentence'].copy()
stcs_diff_MNE  = stcs_MNE['sentence'].copy()
    
stcs_diff_dSPM._data -= stcs_dSPM['scrambled'].data
stcs_diff_MNE._data  -= stcs_MNE['scrambled'].data

#%% ===========================================================================
# DIPOLE FIT - not covered in tutorial
# =============================================================================

evokeds = mne.read_evokeds(join(save_path, 'reading-ave.fif'))
cov = mne.read_cov(join(save_path, 'reading-cov.fif'))
bem = mne.read_bem_solution(join(bem_folder, subject + '-5120-bem-sol.fif'))
trans = join(bem_folder, subject + '-trans.fif')

dipoles = dict()

for evoked in evokeds:
    evoked.crop(0.090, 0.110)
    dipoles[evoked.comment] = dict()
    dipoles[evoked.comment]['dip'], dipoles[evoked.comment]['residual'] = \
        mne.fit_dipole(evoked, cov, bem, trans, n_jobs=1)
        
        
#%% ===========================================================================
# MORPH TO FSAVERAGE   
# =============================================================================

spacing = 'oct6'    
src_name = join(bem_folder, subject + '-' + spacing[:3] + '-' + spacing[-1] + \
                '-src.fif')

src = mne.read_source_spaces(src_name)
stcs_morphs = dict()

for stc in stcs:
    morph = mne.compute_source_morph(src, subject, 'fsaverage', subjects_dir)
    stcs_morphs[stc] = morph.apply(stcs[stc])