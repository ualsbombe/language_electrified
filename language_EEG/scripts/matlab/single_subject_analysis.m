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
                    'language_EEG');
data_dir = fullfile(project_dir, 'data');
figures_dir = fullfile(project_dir, 'figures');
fieldtrip_dir = fullfile(home_dir, 'matlab', 'fieldtrip-20190716');

addpath(fieldtrip_dir);
ft_defaults

%% EXPLORE DATA STRUCTURE

cfg = [];
cfg.dataset = fullfile(data_dir, 'subj2.vhdr');
cfg.trialdef.eventtype = '?';

temp = ft_definetrial(cfg);

%% DEFINE DATA STRUCTURE

orthographic_animals = {'S111' 'S121' 'S131' 'S141'};
orthographic_tools   = {'S151' 'S161' 'S171' 'S181'};

visual_animals       = {'S112' 'S122' 'S132' 'S142'};
visual_tools         = {'S152' 'S162' 'S172' 'S182'};

auditory_animals     = {'S113' 'S123' 'S133' 'S143'};
auditory_tools       = {'S153' 'S163' 'S173' 'S183'};

cfg = [];
cfg.dataset = fullfile(data_dir, 'subj2.vhdr');
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = [orthographic_animals orthographic_tools ...
                           visual_animals visual_tools ...
                           auditory_animals auditory_tools];
cfg.trialdef.prestim    = 0.300; % s
cfg.trialdef.poststim   = 0.700; % s

cfg_def = ft_definetrial(cfg);

%% READ IN AND PREPROCESS DATA

% re-referencing
cfg_def.reref          = 'yes'; % create a new reference
cfg_def.channel        = 'all'; % apply reference to all EEG channels
cfg_def.implicitref    = 'M1'; % the non-recorded reference
cfg_def.refmethod      = 'avg'; % method of referencing
cfg_def.refchannel     = {'M1' '53'}; % the average (see "refmethod" above) of
                                      % these two channels will be the new
                                      % reference; 53 is M2, the recorded
                                      % reference channel

% filtering and demeaning (baselining)
cfg_def.bpfilter       = 'yes'; % apply a band pass filter
cfg_def.bpfreq         = [1 30]; % Hz; range of band pass filter
cfg_def.demean         = 'yes'; % divide the amplitude of each trial by the
                                %
                                % mean amplitude in the time window chosen in
                                % "baselinewindow" below
cfg_def.baselinewindow = [-Inf 0]; % from beginning to time zero

data = ft_preprocessing(cfg_def);

%% EXTRACTING EOG CHANNELS
% create vertical EOG
cfg             = [];
cfg.channel     = {'50' '64'};% electrodes placed vertically around the left eye
cfg.reref       = 'yes';
cfg.implicitref = []; % means no implicit reference
cfg.refchannel  = '50'; % reference '50' and '64' to '50
cfg.refmethod   = 'avg';

eog_vertical = ft_preprocessing(cfg, data);

% keep only one channel and rename it to EEGv
cfg         = [];
cfg.channel = '64'; % the channel to keep

eog_vertical = ft_selectdata(cfg, eog_vertical);
eog_vertical.label = {'EOGv'}; % renaming it

% create horizontal EOG (logic similar to above)
cfg             = [];
cfg.channel     = {'51' '60'};
cfg.reref       = 'yes';
cfg.implicitref = [];
cfg.refchannel  = '51';
cfg.refmethod   = 'avg';

eog_horizontal = ft_preprocessing(cfg, data);

% keep only channel and rename it to EEGh (logic similar to above)
cfg         = [];
cfg.channel = '60';

eog_horizontal = ft_selectdata(cfg, eog_horizontal);
eog_horizontal.label = {'EOGh'};

% rename '53' to M2 for consistency
channel_1_index = strcmp(data.label, '53');
data.label{channel_1_index} = 'M2';

% keep only the EEG channels from the original data
cfg         = [];
cfg.channel = setdiff(1:60, [50 51 60 64]); % removing channels 50, 51, 60 & 64
                                            % since they do not contain EEG
                                            % signal, but have been
                                            % converted to EOG
data = ft_selectdata(cfg, data);

% finally append EOGv and EOGh data to EEG data
cfg = [];
data = ft_appenddata(cfg, data, eog_vertical, eog_horizontal);

%% REJECT TRIALS BASED ON OBJECTIVE THRESHOLD (150 µV)

threshold = 150; % µV
n_trials = length(data.trial);
bad_indices = []; % initialize empty - we don't know how many bad trials
                        % there will be
                                              
for trial_index = 1:n_trials
    
    trial = data.trial{trial_index};  % get trial n
    % check whether threshold is exceeed in either the negative or (|) the
    % positive direction for any of the time points on any of the
    % electrodes
    threshold_exceeded = any(any(trial > threshold | trial < -threshold));
    if threshold_exceeded
        % add the index of each trial that exceeded the threshold
        bad_indices = [bad_indices trial_index]; %#ok<AGROW>
    end
    
end

% keep only the good indices
good_indices = setdiff(1:n_trials, bad_indices);

cfg =        [];
cfg.trials = good_indices;

cleaned_data = ft_selectdata(cfg, data);

% and make common reference

cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';

cleaned_data = ft_preprocessing(cfg, cleaned_data);
    
%% EVENT-RELATED RESPONSES

categories = {
                [111 121 131 141] % orthographic_animals 
                [151 161 171 181] % orthographic_tools   

                [112 122 132 142] % visual_animals       
                [152 162 172 182] % visual_tools 

                [113 123 133 143] % auditory_animals
                [153 163 173 183] % auditory_tools
             };

n_categories = length(categories);
ERPs = cell(1, n_categories);

for category_index = 1:n_categories
    
    category = categories{category_index};
    
    cfg                  = [];
    cfg.trials           = ismember(cleaned_data.trialinfo, category);
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'prestim';
   
    ERPs{category_index} = ft_timelockanalysis(cfg, cleaned_data);
    
end

%% DIFFERENCE WAVES

categories = {
                [111 121 131 141] % orthographic_animals 
                [151 161 171 181] % orthographic_tools   

                [112 122 132 142] % visual_animals       
                [152 162 172 182] % visual_tools 

                [113 123 133 143] % auditory_animals
                [153 163 173 183] % auditory_tools
             };

n_categories = length(categories);
differences_semantic_categories = cell(1, n_categories/2);

for category_index = 2:2:n_categories
    
    cfg = [];
    cfg.operation = 'x1 - x2'; % subtract animals from tools
    cfg.parameter = 'avg';
    
    differences_semantic_categories{category_index/2} = ...
        ft_math(cfg, ERPs{category_index}, ERPs{category_index-1});
    
end

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 35, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)
        
%% PLOT ERPs
% multiplot
mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg             = [];
cfg.layout      = 'easycapM10.mat';
cfg.xlim        = [-0.300 0.698]; % s
cfg.ylim        = [-10 10]; % µV
cfg.limittext   = sprintf(['Auditory Response\nx: Time: ' ...
                        num2str(cfg.xlim(1), '%0.3f') '-' ...
                        num2str(cfg.xlim(2), '%0.3f') ...
                        ' s\n\t\ty: Electric Potential ' ...
                        num2str(cfg.ylim(1), '%0.1f') '-' ...
                        num2str(cfg.ylim(2), '%0.1f') ' µV']);
cfg.fontsize    = 25;
cfg.showscale   = 'no';

ft_multiplotER(cfg, ERPs{5});

%single channel plot
splot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg           = [];
cfg.xlim      = [-0.300 0.698]; % s
cfg.ylim      = [-10 10]; % µV
cfg.layout    = 'easycapM10.mat';
cfg.channel   = '1';
cfg.title     = sprintf(['Auditory Response\nCentral Electrode (' ...
                         cfg.channel ')']);
cfg.fontsize  = 35;
cfg.linewidth = 3;
xlabel('Time (s)')
ylabel('Electric Potential (µV)')

ft_singleplotER(cfg, ERPs{5});
hold on
plot([cfg.xlim(1) cfg.xlim(2)], [0 0], 'k--'); % horizontal line
plot([0 0], [cfg.ylim(1) cfg.ylim(2)], 'k--'); % vertical line

% topoplot
tplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg          = [];
cfg.layout   = 'easycapM10.mat';
cfg.xlim     = [0.120 0.120]; % s
cfg.zlim     = [-4 4]; % µV
cfg.comment  = sprintf(['Auditory N1\nTime: ' ...
                        num2str(cfg.xlim(1), '%0.3f') ' s']);
cfg.fontsize = 25;
c = colorbar;
c.Label.String = 'Electric Potential (µV)';
c.Label.FontSize = 40;

ft_topoplotER(cfg, ERPs{5});

%% SAVE PLOTS
print(mplot, fullfile(figures_dir, 'multiplot_auditory.jpeg'), '-djpeg', ...
    '-r300')
print(splot, fullfile(figures_dir, 'singleplot_auditory.jpeg'), '-djpeg', ...
    '-r300')
print(tplot, fullfile(figures_dir, 'topoplot_auditory.jpeg'), '-djpeg', ...
    '-r300')

%% PLOT ERPs (ROUGH)
% multiplot
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg             = [];
cfg.layout      = 'easycapM10.mat';

ft_multiplotER(cfg, ERPs{5});
%single channel plot
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg           = [];

cfg.layout    = 'easycapM10.mat';
cfg.channel   = '1';

ft_singleplotER(cfg, ERPs{5});

% topoplot
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
cfg          = [];
cfg.layout   = 'easycapM10.mat';
cfg.xlim     = [0.120 0.120]; % s

ft_topoplotER(cfg, ERPs{5});

%% PLOT DIFFERENCE WAVES

cfg        = [];
cfg.layout = 'easycapM10.mat';

figure;
ft_multiplotER(cfg, differences_semantic_categories{2});

%% SOURCE RECONSTRUCTION (load templates)

headmodel = ft_read_headmodel('standard_bem.mat');
elec = ft_read_sens('easycap-M10.txt');
sourcemodel = ft_read_headshape('cortex_5124.surf.gii');

% realign electrodes
cfg           = [];
cfg.method    = 'project';
cfg.headshape = headmodel.bnd(1); % head surface

elec_realigned = ft_electroderealign(cfg, elec);

% get everything in SI units

headmodel = ft_convert_units(headmodel, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');
elec = ft_convert_units(elec, 'm');
elec_realigned = ft_convert_units(elec_realigned, 'm');

% plot models

close all

quality_plot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(1, 2, 1)
hold on
ft_plot_headmodel(headmodel, 'facealpha', 0.3)
ft_plot_sens(elec, 'facecolor', 'red', 'label', 'label', ...
    'fontsize', 30, 'elecsize', 40)
ft_plot_mesh(sourcemodel)
view(-90, 0)
title('No projection');

subplot(1, 2, 2)
hold on
ft_plot_headmodel(headmodel, 'facealpha', 0.3)
ft_plot_sens(elec_realigned, 'facecolor', 'red', 'label', 'label', ...
    'fontsize', 30, 'elecsize', 40)
ft_plot_mesh(sourcemodel)
view(-90, 0)
title('Projected onto scalp')

%% SAVE QUALITY PLOT

print(quality_plot, fullfile(figures_dir, 'quality_plot.jpeg'), '-djpeg', ...
    '-r300')

%% LEADFIELD

cfg = [];
cfg.elec = elec_realigned;
cfg.sourcemodel.pos = sourcemodel.pos;
cfg.sourcemodel.inisde = 1:size(sourcemodel.pos, 1);
cfg.headmodel = headmodel;

leadfield = ft_prepare_leadfield(cfg);

%% PLOT LEADFIELD

close all

channel_1 = '1';
channel_2 = '43';
% find channel indices
channel_1_index = find(strcmp(elec_realigned.label, channel_1));
channel_2_index = find(strcmp(elec_realigned.label, channel_2));

pos_1 = elec_realigned.chanpos(channel_1_index, :); %% ... and their positions
pos_2 = elec_realigned.chanpos(channel_2_index, :);

% find the distances between the sensor and all sources
distances_1 = zeros(1, length(leadfield.inside)); 
distances_2 = distances_1;

% in the loop we go through the sources one-by-one ...
for index = 1:length(distances_1)     
    distances_1(index) = norm(pos_1 - leadfield.pos(index, :));
    distances_2(index) = norm(pos_2 - leadfield.pos(index, :));
end

% only take the ones that are actually inside 
distances_1 = distances_1(leadfield.inside);
distances_2 = distances_2(leadfield.inside);

[~, source_index_1] = min(distances_1);
[~, source_index_2] = min(distances_2);

inside_indices = find(leadfield.inside);

source_indices = [source_index_1 source_index_2];

titles = {'Central Source' 'Posterior source'};

leadfield_plot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

moment = 10e-9; % 10 nAm
scaling = 1e6; % to make sure ft_plot_sens can plot the colour differences

for ii = 1:length(source_indices)
    
    subplot(1, 2, ii)

    source_index = inside_indices(source_indices(ii));
    % find the leadfield corresponding to the source closest to the sensor
    this_leadfield = leadfield.leadfield{source_index};  

    n_sensors = length(this_leadfield);
    channels_leadfield = zeros(1, n_sensors);

    % getting the field strength for each leadfield vector
    for n_sensor = 1:n_sensors
        channels_leadfield(n_sensor) = norm(this_leadfield(n_sensor, :));
    end
    
    % preparing for colouring the sensors with different
    % hues of red according to field strength
    colours = ones(n_sensors, 3);
    for colour_index = 1:length(colours)
            % 1e4 is chosen to make the value sufficiently bigger than 0,
            % such that colour changes can be seen
            colours(colour_index, :) = ...
                [channels_leadfield(colour_index) * 1e0 0 0] * moment * scaling; 
    end

    ft_plot_sens(elec_realigned, 'facealpha', 0.8, 'facecolor', colours, ...
        'elecsize', 80); %% plot the sensors with the colouring above
    hold on
%     ft_plot_topo3d(elec.chanpos, channels_leadfield, 'facealpha', 0.8);

    source_colours = zeros(length(leadfield.inside), 3);
    source_colours(source_index, :) = [0 0 1]; %% plot the nearest source blue
    vertexsizes = 10 * ones(length(leadfield.inside));
    vertexsizes(source_index) = 250; % visualize it as a bigger than the others

    ft_plot_mesh(leadfield, 'vertexcolor', source_colours, ...
        'facealpha', 0.8, 'vertexsize', vertexsizes)
    ft_plot_headmodel(headmodel, 'facealpha', 0.8);
    view(-90, 0)
    title(['Lead field: ' titles{ii}])
end

colour_matrix = zeros(101, 3);
colour_matrix(:, 1) = 0:0.01:1;
colormap(colour_matrix);

c = colorbar;
colormap(colour_matrix);
c.Label.String = 'Sensitivity (low-to-high, a.u.)'; 


%% PRINT LEADFIELD PICTURE

print(leadfield_plot, fullfile(figures_dir, 'leadfield_plot.jpeg'), ...
    '-djpeg', '-r300')

%% INVERSE SOLUTIONS 

MNEs = cell(1, n_categories);

for category_index = 1:n_categories
    
    cfg                    = [];
    cfg.method             = 'mne';
    cfg.senstype           = 'eeg';
    cfg.elec               = elec_realigned;
    cfg.channel            = cleaned_data.label(1:57);
    cfg.sourcemodel        = leadfield;
    cfg.headmodel          = headmodel;
    cfg.mne.prewhiten      = 'no';
    cfg.mne.lambda         = 0.1;
    cfg.mne.scalesourcecov = 'yes';
    
    MNE = ft_sourceanalysis(cfg, ERPs{category_index});
    MNE.tri = sourcemodel.tri;
    
    MNEs{category_index} = MNE;
    
end

%% LAMBDA TEST

lambdas = [0.1 0.5 1 2 3 10 50 100];
% n_lambdas = length(lambdas);
% 
% MNEs = cell(1, n_lambdas);
% 
% for lambda_index = 1:n_lambdas
% 
%     cfg                    = [];
%     cfg.method             = 'mne';
%     cfg.senstype           = 'eeg';
%     cfg.elec               = elec_realigned;
%     cfg.channel            = cleaned_data.label(1:57);
%     cfg.sourcemodel        = leadfield;
%     cfg.headmodel          = headmodel;
%     cfg.mne.prewhiten      = 'no';
%     cfg.mne.lambda         = lambdas(lambda_index);
%     cfg.mne.scalesourcecov = 'yes';
%     
%     MNE = ft_sourceanalysis(cfg, ERPs{5});
%     MNE.tri = sourcemodel.tri;
%     
%     MNEs{lambda_index} = MNE;
% 
% end
%% VISUAL PLOT

% for lambda_index = 1:n_lambdas
% 
%     cfg = [];
%     cfg.method = 'surface';
%     cfg.funparameter = 'pow';
%     cfg.funcolorlim = [0 3e-9]; % Am
%     cfg.funcolormap = 'parula';
%     cfg.latency = 0.120; % s
% 
%     ft_sourceplot(cfg, MNEs{lambda_index});
% end

%% MNE PLOT

close all

category_index = 5;

this_mne = MNEs{category_index};
this_mne.avg.pow_nAm = this_mne.avg.pow * 1e9; 

cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow_nAm';
cfg.funcolorlim = [0 50]; % nAm
cfg.funcolormap = 'parula';
cfg.latency = 0.100; % s
cfg.colorbar = 'no';

ft_sourceplot(cfg, this_mne);
view(90, 0)

c = colorbar;
c.Label.String = 'Current (nAm)';
mne_plot = gcf;
set(mne_plot, 'units', 'normalized', 'outerposition', [0 0 1 1]);
title(['N100 localization (' num2str(1e3 * cfg.latency) ' ms)']);

%% SAVE MNE PLOT

print(mne_plot, fullfile(figures_dir, 'mne_plot.jpeg'), ...
    '-djpeg', '-r300')

%% LEADFIELD GRID

cfg = [];
cfg.elec = elec_realigned;
cfg.headmodel = headmodel;
cfg.resolution = 0.01;
cfg.sourcemodel.unit = 'm';
cfg.channel = cleaned_data.label(1:57);

sourcemodel = ft_prepare_leadfield(cfg);

%% DIPOLE FITTING

cfg = [];
cfg.latency = [0.110 0.130];
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel;
cfg.model = 'regional';
cfg.nonlinear = 'yes';
cfg.gridsearch = 'yes';
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.elec = elec_realigned;

dipole_fit = ft_dipolefitting(cfg, ERPs{5});

%% PLOT DIPOLE FIT

load('standard_mri.mat')
mri = ft_convert_units(mri, 'm');
dipole_mri_plot = figure;
hold on
ft_plot_dipole(dipole_fit.dip.pos(1, :), ...
                mean(dipole_fit.dip.mom(1:3, :), 2), 'unit', 'm');
ft_plot_dipole(dipole_fit.dip.pos(2, :), ...
                mean(dipole_fit.dip.mom(4:6, :), 2), 'unit', 'm');

pos = mean(dipole_fit.dip.pos, 1);

orientations = {[1 0 0] [0 1 0] [0 0 1]};
n_orientations = length(orientations);

for orientation_index = 1:n_orientations
    
    orientation = orientations{orientation_index};
    
    ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', pos, ...
                  'orientation', orientation);
end

axis tight
axis off

%% SAVE DIPOLE MRI PLOT

print(dipole_mri_plot, fullfile(figures_dir, 'dipole_mri_plot.jpeg'), ...
    '-djpeg', '-r300')

%% REVERT TO SI UNITS

cfg = [];
cfg.operation = 'x1 / 1e6';
cfg.parameter = 'avg';

ERP_SI = ft_math(cfg, ERPs{5});

%% ESTIMATE DIPOLE TIME COURSES

cfg = [];
cfg.latency = 'all';
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel;
cfg.model = 'regional';
cfg.nonlinear = 'yes';
cfg.gridsearch = 'no';
cfg.numdipoles = 2;
cfg.symmetry = 'x';
cfg.elec = elec_realigned;
cfg.dip.pos = dipole_fit.dip.pos;

dipole_time_courses = ft_dipolefitting(cfg, ERPs{5});

%% PLOT ESTIMATED TIME COURSES

dipole_tc_plot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

n_times = length(dipole_time_courses.time);
amplitudes = zeros(2, n_times);

for time_index = 1:n_times
    amplitudes(1, time_index) = ...
                            norm(dipole_time_courses.dip.mom(1:3, time_index));
    amplitudes(2, time_index) = ...
                            norm(dipole_time_courses.dip.mom(4:6, time_index));
end

amplitudes = amplitudes / 1e6; %% Put it in SI units

plot(dipole_time_courses.time, amplitudes);
xlim([dipole_time_courses.time(1), dipole_time_courses.time(end)]);
xlabel('Time (s)')
ylabel('Amplitude (Am)')
title('Dipole time courses')
legend({'Right Hemisphere', 'Left Hemisphere'})

%% SAVE DIPOLE TIME COURSES

print(dipole_tc_plot, fullfile(figures_dir, 'dipole_tc_plot.jpeg'), ...
    '-djpeg', '-r300');