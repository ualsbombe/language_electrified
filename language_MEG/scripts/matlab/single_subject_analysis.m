%% CLEAN PATH AND WORKSPACE

clear variables
restoredefaultpath

% change "home_dir", "data_dir", project_dir and "fieldtrip_dir" such that
% it applies to your setup.
% If you prefer having the "fieldtrip_dir" on your default path
% delete the lines about "fieldtrip_dir" and adding "fieldtrip_dir" to your
% path. Do remember to run "ft_defaults". This makes sure that the
% appropriate FieldTrip subdirectories are loaded
home_dir = '/home/lau/';
project_dir = fullfile(home_dir, 'analyses', 'fieldtrip_datasets', ...
                    'language_MEG');
data_dir = fullfile(project_dir, 'data');
figures_dir = fullfile(project_dir, 'figures');
fieldtrip_dir = fullfile(home_dir, 'matlab', 'fieldtrip-20190716');
derived_meg_dir = fullfile(data_dir, 'sub-V1002', 'meg', 'derived');
derived_mr_dir = fullfile(data_dir, 'sub-V1002', 'anat', 'derived');

addpath(fieldtrip_dir);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 35, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)

%% DEFINE DATA STRUCTURE

sentence  = {2 4};
scrambled = {6 8};

cfg = [];
cfg.dataset = fullfile(data_dir, 'sub-V1002', 'meg', ...
    'sub-V1002_task-visual_meg.ds');
cfg.trialdef.eventtype  = 'frontpanel trigger';
cfg.trialdef.eventvalue = [sentence scrambled];
cfg.trialdef.prestim    = 0.514; % s
cfg.trialdef.poststim   = 1.036; % s % allowing for adjsutment

cfg_def = ft_definetrial(cfg);

%% READ IN DATA AND ADJUST TIMELINE

cfg_def.continuous = 'yes'; %% the data was recorded continuously

data = ft_preprocessing(cfg_def);

cfg = [];
cfg.offset = -0.036 * data.fsample; %% push forward 36 ms in time

data = ft_redefinetrial(cfg, data);

save(fullfile(derived_meg_dir, 'preprocessed_data'),  '-struct', 'data');

%% GET TRIAL INDICES

n_events = 2;
trials = cell(1, n_events);

for event_index = 1:n_events
    
    if event_index == 1
        trial_values = sentence;
    else
        trial_values = scrambled;
    end
    
    n_trials = length(data.trialinfo);
    trial_indices = zeros(1, n_trials);
    
    for trial_index = 1:n_trials
        if any(data.trialinfo(trial_index) == cell2mat(trial_values))
            trial_indices(trial_index) = 1; 
        end
    end
   
    trials{event_index} = logical(trial_indices); %% make boolean and add
                                                  % to cell array
    
end
    
%% EVENT-RELATED FIELDS (ERFs)

n_events = 2;
ERFs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    cfg.preproc.demean = 'yes'; %% demean the data
    cfg.preproc.baselinewindow = [-0.150 0]; %% use this time window to demean
    cfg.preproc.lpfilter = 'yes'; %% apply low-pass filter ...
    cfg.preproc.lpfreq = 40; %% at this frequency
    cfg.trials = trials{event_index}; %% use only these trials
    cfg.latency = [-0.150 0.500]; %% in only this time window
    cfg.covariance = 'yes';
    cfg.covariancewindow = cfg.preproc.baselinewindow;
    cfg.channel = 'MEG';
    
    ERFs{event_index} = ft_timelockanalysis(cfg, data);
end

save(fullfile(derived_meg_dir, 'ERFs'), 'ERFs');

%% PLOT TIMELOCKEDS

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim = [0.100 0.100]; %% plot the visual response at 100 ms
cfg.zlim = [-120e-15 120e-15];
cfg.comment = ' ';

mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(1, 2, 1)
ft_topoplotER(cfg, ERFs{1});
title('Sentence (100 ms)')

subplot(1, 2, 2)
ft_topoplotER(cfg, ERFs{2});
title('Scrambled (100 ms)')
c = colorbar;
c.Label.String = 'Magnetic Field (T)';

print(mplot, fullfile(figures_dir, 'topoplot_fieldtrip.jpeg'), '-djpeg', ...
    '-r300')

%% BASELINE DATA

cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-Inf Inf];

baselined_data = ft_preprocessing(cfg, data);

%% DETAILED TFRs

n_events = 2;
detailed_tfrs = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'mtmconvol'; %% run tfr
    cfg.taper = 'hanning'; %% use a hanning taper
    cfg.output = 'pow'; %% return power
    cfg.foi = 2:2:40; %% frequencies to centre on
    cfg.t_ftimwin = ones(1, length(cfg.foi)) * 0.500; %% use a sliding window
                                                      % of 400 ms 
    cfg.toi = -0.100:(1/baselined_data.fsample):0.500; %% time points to
                                                       % centre on
    cfg.pad = 'nextpow2';
    cfg.channel = 'meg'; %% which channels to use
    cfg.trials = trials{event_index}; %% which condition, sentence or scrambled
    
    detailed_tfrs{event_index} = ft_freqanalysis(cfg, baselined_data);
end

save(fullfile(derived_meg_dir, 'detailed_tfrs'), 'detailed_tfrs');

%% PLOT DETAILED TFRs

n_events = 2;

cfg = [];
cfg.baseline = [-Inf 0]; %% baseline window
cfg.baselinetype = 'relative'; %% for each frequency demean based on the window
                                % above
cfg.layout = 'CTF275.lay';
cfg.zlim = [0.7 1.6];
cfg.channel = 'MLT15';
cfg.colorbartext = 'Power relative to baseline';
cfg.fontsize = 35; %% of title

names = {'Sentence' 'Scrambled'};
mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for event_index = 1:n_events
    
    if event_index == 1
        cfg.colorbar = 'no';
    elseif event_index == 2
        cfg.colorbar = 'yes';
    end
    subplot(1, 2, event_index);
    cfg.title = names{event_index};
    ft_singleplotTFR(cfg, detailed_tfrs{event_index});
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    
end

print(mplot, fullfile(figures_dir, 'detailed_tfr.jpeg'), '-djpeg', '-r300')

%% PLOT DETAILED TFRs (no baseline)

n_events = 2;

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.channel = 'MLT15';
cfg.colorbartext = 'Power (T^2)';
cfg.fontsize = 35; %% of title

names = {'Sentence' 'Scrambled'};
mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for event_index = 1:n_events
    
    if event_index == 1
        cfg.colorbar = 'no';
    elseif event_index == 2
        cfg.colorbar = 'yes';
    end
    subplot(1, 2, event_index);
    cfg.title = names{event_index};
    ft_singleplotTFR(cfg, detailed_tfrs{event_index});
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    
end

print(mplot, fullfile(figures_dir, 'detailed_tfr_no_baseline.jpeg'), ...
    '-djpeg', '-r300')

%% MR PREPROCESSING

% read in
mri = ft_read_mri(fullfile(data_dir, 'sub-V1002', 'anat', ...
                           'sub-V1002_space-CTF_T1w.nii'));

mri.coordsys = 'ctf';

cfg = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri);
                     
cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

segmented_mri = ft_volumesegment(cfg, mri);

save(fullfile(derived_mr_dir, 'mri_resliced'), 'mri_resliced');
save(fullfile(derived_mr_dir, 'segmented_mri'), 'segmented_mri')

%% CONSTRUCT MESH AND HEADMODEL

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, segmented_mri);

cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, segmented_mri);

% construct headmodel

cfg = [];
cfg.method = 'singleshell';

headmodel = ft_prepare_headmodel(cfg, mesh_brain);
headmodel = ft_convert_units(headmodel, 'm');

save(fullfile(derived_mr_dir, 'headmodel'), 'headmodel')

%% PLOT ALIGNMENT (QUALITY CHECK)

sens = ft_read_sens(fullfile(data_dir, 'sub-V1002', 'meg', ...
    'sub-V1002_task-visual_meg.ds'), 'senstype', 'meg');
headshape = ft_read_headshape(fullfile(data_dir, 'sub-V1002', 'meg', ...
    'sub-V1002_headshape.pos'));
mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
ft_plot_sens(sens, 'unit', 'm');
ft_plot_headshape(headshape, 'unit', 'm', 'vertexcolor', 'red', ...
    'vertexsize', 50);
ft_plot_mesh(mesh_scalp, 'unit', 'm', 'facealpha', 0.8);
ft_plot_axes([], 'unit', 'm')

view(180, 0)
zoom(1.5)

print(mplot, fullfile(figures_dir, 'quality_check.jpeg'), '-djpeg', ...
    '-r300')

%% LEADFIELD

sens = ft_read_sens(fullfile(data_dir, 'sub-V1002', 'meg', ...
    'sub-V1002_task-visual_meg.ds'), 'senstype', 'meg');
sens = ft_convert_units(sens, 'm');

cfg = [];
cfg.grad = sens;
cfg.headmodel = headmodel;
cfg.resolution = 0.01;
cfg.sourcemodel.unit = 'm';
cfg.channel = 'MEG';

sourcemodel = ft_prepare_leadfield(cfg);

save(fullfile(derived_mr_dir, 'sourcemodel'), 'sourcemodel')

%% PLOT SOURCE MODEL

mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
inside_pos = sourcemodel.pos(sourcemodel.inside, :);
outside_pos = sourcemodel.pos(~sourcemodel.inside, :);
ft_plot_headmodel(headmodel, 'facealpha', 0.6)
ft_plot_mesh(inside_pos, 'vertexcolor', 'blue', 'vertexsize', 30)
ft_plot_mesh(outside_pos, 'vertexcolor', 'red', 'vertexsize', 30);

view(33, 55)

print(mplot, fullfile(figures_dir, 'sourcemodel.jpeg'), '-djpeg', '-r300');

%% TFR FOR BEAMFORMER

n_events = 2;
tfrs = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'mtmconvol'; %% run tfr
    cfg.taper = 'hanning'; %% use a hanning taper
    cfg.output = 'powandcsd'; %% return both power and cross-spectral density
    cfg.t_ftimwin = [0.400 0.400]; %% use a sliding window of 200 ms
    cfg.toi = -0.100:0.050:0.500; %% time points to centre on
    cfg.pad = 'nextpow2'; %% which padding to use
    cfg.channel = 'meg'; %% which channels to use
    cfg.trials = trials{event_index}; %% which condition, sentence or scrambled
    cfg.foi = [5 10]; %% the two frequencies focused on

    tfrs{event_index} = ft_freqanalysis(cfg, baselined_data);
end

cfg = rmfield(cfg, 'trials');
tfr_combined = ft_freqanalysis(cfg, baselined_data);

tfrs{3} = tfr_combined;

for tfr_index = 1:length(tfrs)
    tfrs{tfr_index}.grad = ft_convert_units(tfrs{tfr_index}.grad, 'm');
end

save(fullfile(derived_meg_dir, 'tfrs'), 'tfrs');

%% PLOT TIME-FREQUENCY ANALYSIS COARSE

n_events = 2;

cfg = [];
cfg.baseline = [-Inf 0]; %% baseline window
cfg.baselinetype = 'relative'; %% for each frequency demean based on the window
                                % above
cfg.layout = 'CTF275.lay';
cfg.zlim = [0.7 1.6];
cfg.channel = 'MLT15';
cfg.colorbartext = 'Power relative to baseline';
cfg.fontsize = 35; %% of title

names = {'Sentence' 'Scrambled'};
mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for event_index = 1:n_events
    
    if event_index == 1
        cfg.colorbar = 'no';
    elseif event_index == 2
        cfg.colorbar = 'yes';
    end
    subplot(1, 2, event_index);
    cfg.title = names{event_index};
    ft_singleplotTFR(cfg, tfrs{event_index});
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    
end

print(mplot, fullfile(figures_dir, 'coarse_tfr.jpeg'), '-djpeg', '-r300')

%% LOAD FOR BEAMFORMER

load(fullfile(derived_meg_dir, 'tfrs.mat'))
load(fullfile(derived_mr_dir, 'sourcemodel.mat'))
load(fullfile(derived_mr_dir, 'headmodel.mat'))

%% BEAMFORMER THETA

n_events = 2;
n_times = length(tfrs{1}.time);
beamformers = cell(1, n_events);

for time_index = 1:n_times
    
    time = tfrs{1}.time(time_index);
    
    cfg = [];
    cfg.method = 'dics'; % use Dynamic Imaging of Coherent Sources
    cfg.frequency = tfrs{1}.freq(1); % Choose the theta band
    cfg.sourcemodel = sourcemodel; 
    cfg.headmodel = headmodel;
    cfg.dics.projectnoise = 'yes'; % project out noise
    cfg.dics.keepfilter = 'yes'; % keep the created spatial filter
    cfg.dics.realfilter = 'yes'; % use the real part of the filter
    cfg.latency = time; % the time 
    
    beamformer_combined = ft_sourceanalysis(cfg, tfrs{3}); % do the combined
    
    cfg.sourcemodel.filter = beamformer_combined.avg.filter; % add filter to
                                                             % sourcemodel
    
    % do for sentence and scrambled
    for event_index = 1:n_events
        
        if time_index == 1
            beamformers{event_index} = cell(1, n_times);
        end
        
        beamformers{event_index}{time_index} = ft_sourceanalysis(cfg, ...
                                                            tfrs{event_index});
                                                        
    end
    
end

save(fullfile(derived_meg_dir, 'beamformers'), 'beamformers');

%% BEAMFORMER THETA (lambda)

n_events = 2;
n_times = length(tfrs{1}.time);
beamformers_lambda = cell(1, n_events);

for time_index = 1:n_times
    
    time = tfrs{1}.time(time_index);
    
    cfg = [];
    cfg.method = 'dics';
    cfg.frequency = tfrs{1}.freq(1);
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = headmodel;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.keepfilter = 'yes';
    cfg.dics.realfilter = 'yes';
    cfg.dics.lambda = '10%';
    cfg.latency = time;
    
    beamformer_combined = ft_sourceanalysis(cfg, tfrs{3});
    
    cfg.sourcemodel.filter = beamformer_combined.avg.filter;
    
    for event_index = 1:n_events
        
        if time_index == 1
            beamformers_lambda{event_index} = cell(1, n_times);
        end
        
        beamformers_lambda{event_index}{time_index} = ft_sourceanalysis(cfg, ...
                                                            tfrs{event_index});
                                                        
    end
    
end

save(fullfile(derived_meg_dir, 'beamformers_lambda'), 'beamformers_lambda');

%% LOAD BEAMFORMERS

load(fullfile(derived_meg_dir, 'beamformers.mat'))
load(fullfile(derived_meg_dir, 'beamformers_lambda.mat'))
load(fullfile(derived_mr_dir, 'mri_resliced.mat'))

%% SET SOURCE PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 20, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)

%% INTERPOLATE WITHOUT CONTRAST AND PLOT

cfg = [];
cfg.parameter = 'pow';
cfg.downsample = 2;

% {1} is sentence, {8} is 250 ms
beamformer_interpolated = ft_sourceinterpolate(cfg, beamformers{1}{8}, ...
                                               mri_resliced);
% SOURCE PLOT

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.crosshair = 'no';
cfg.comment = ' ';
cfg.colorbartext = 'Source Power (Am)^2';

ft_sourceplot(cfg, beamformer_interpolated);

mplot = gcf;
set(mplot, 'units', 'normalized', 'outerposition', [0 0 1 1])

print(fullfile(figures_dir, 'beamformer.jpg'), '-djpeg', '-r300');
                                           
%% CONTRAST, INTERPOLATE AND PLOT

% {1} is sentence, {2} is scrambled, {8} is 250 ms

% how much greater/smaller is the power difference between the two
% conditions to the power in the scrambled condition
cfg = [];
cfg.operation = '(x1 - x2) / x2';
cfg.parameter = 'pow';
                                                                            
beamformer_contrast = ft_math(cfg, beamformers{1}{8}, beamformers{2}{8});   
                                                                                           
cfg = [];
cfg.parameter = 'pow';
cfg.downsample = 2;

beamformer_contrast_interpolated = ft_sourceinterpolate(cfg, ...
                                                        beamformer_contrast, ...
                                                        mri_resliced);
                                                    
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.crosshair = 'yes';
cfg.comment = ' ';
cfg.colorbartext = 'Source Power Ratio';
cfg.funcolorlim = [0.10 0.3];
cfg.location = [28 53 63];

ft_sourceplot(cfg, beamformer_contrast_interpolated);

mplot = gcf;
set(mplot, 'units', 'normalized', 'outerposition', [0 0 1 1])

print(fullfile(figures_dir, 'beamformer_contrast.jpg'), '-djpeg', '-r300'); 

%% CONTRAST, INTERPOLATE AND PLOT (lambda)

% {1} is sentence, {2} is scrambled, {8} is 250 ms

% how much greater/smaller is the power difference between the two
% conditions to the power in the scrambled condition
cfg = [];
cfg.operation = '(x1 - x2) / x2';
cfg.parameter = 'pow';
                                                                            
beamformer_contrast = ft_math(cfg, beamformers_lambda{1}{8}, ...
                                beamformers_lambda{2}{8});   
                                                                                           
cfg = [];
cfg.parameter = 'pow';
cfg.downsample = 2;

beamformer_contrast_interpolated = ft_sourceinterpolate(cfg, ...
                                                        beamformer_contrast, ...
                                                        mri_resliced);
                                                    
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.crosshair = 'yes';
cfg.comment = ' ';
cfg.colorbartext = 'Source Power Ratio';
cfg.funcolorlim = [0.1 0.6];
cfg.location = [28 53 63];

ft_sourceplot(cfg, beamformer_contrast_interpolated);

mplot = gcf;
set(mplot, 'units', 'normalized', 'outerposition', [0 0 1 1])

print(fullfile(figures_dir, 'beamformer_contrast_lambda.jpg'), '-djpeg', ...
    '-r300'); 

%% WARPED SOURCE MODEL

load(fullfile(derived_mr_dir, 'headmodel.mat'))

sens = ft_read_sens(fullfile(data_dir, 'sub-V1002', 'meg', ...
    'sub-V1002_task-visual_meg.ds'), 'senstype', 'meg');
sens = ft_convert_units(sens, 'm');

mri = ft_read_mri(fullfile(data_dir, 'sub-V1002', 'anat', ...
                           'sub-V1002_space-CTF_T1w.nii'));

mri.coordsys = 'ctf';

MNI_template = fullfile(fieldtrip_dir, 'template', 'sourcemodel', ...
                        'standard_sourcemodel3d10mm.mat');                 

cfg = [];
cfg.grad = sens;
cfg.headmodel = headmodel;
cfg.sourcemodel.unit = 'm';
cfg.sourcemodel.warpmni = 'yes';
cfg.sourcemodel.template = MNI_template;
cfg.sourcemodel.nonlinear = 'yes';
cfg.mri = mri;
cfg.channel = 'MEG';

warped_sourcemodel = ft_prepare_leadfield(cfg);

mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on
ft_plot_headmodel(headmodel)
ft_plot_mesh(warped_sourcemodel)
view(180, 0)

print(mplot, fullfile(figures_dir, 'warped_sourcemodel.jpeg'), '-djpeg', ...
    '-r300')