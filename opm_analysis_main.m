% OPM Pipeline adapted from Christophs opm_general
%
% Make sure mri_main.m is run prior to running this script!

%% Reset all
clear all
close all
restoredefaultpath
addpath("C:\Github\WM_Distractors_OPM\Matlab")
A_Initialize_environment_OPM
global STATIC;
STATIC = StaticVarClass_WM_shift_OPM();


%% Base paths
server = false;
% Laptop:
base_data_path = 'G:\SV10_OPM_Distractor_Processing\SV10_RawData\MEG\raw';
base_save_path = 'G:\SV10_OPM_Distractor_Processing\ProcessedData_Natmeg';
base_matlab_path = 'C:\Github\WM_Distractors_OPM\Matlab';
project_scripts_path = 'C:\toolbox\opm_general';

%% Set up fieldtrip
% addpath(fullfile(base_matlab_path,'fieldtrip')) % Fieldtrip path
% addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Analyse squid data?
squid = false

%% Overwrite
overwrite = [];
overwrite.preproc = true;
overwrite.timelock = true;
overwrite.coreg = true;
overwrite.dipole = false;
overwrite.mne = false;


%% Params
params = [];
params.pre = 1.5; % Trial prestim in seconds
params.post = 9.; % Trial poststim in seconds
params.pad = 0.2; % Trial (pre and post) padding in seconds
params.delay = 0.01; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

params.filter = [];
params.filter.hp_freq = [1]; % Highpass cutoff frequency
params.filter.lp_freq = [150]; % Lowpass cutoff frequency
%params.filter.bp_freq = [1 50]; % Bandpass cutoff frequencies
params.filter.notch = sort([50 60]); % Notch (bandstop) filter frequencies

params.ds_freq = 250; % Downsample frequency. If empty or not defined no downsampling will be applied

% Spatiotemporal filter (OPM-MEG only)
params.do_hfc = false;
params.hfc_order = 2;
params.do_amm = true;
params.amm_in = 12;
params.amm_out = 3;
params.amm_thr = 0.99;

params.n_comp = 40; % Number of ICA components
params.manual_ica = false; % Manually select ICA components to remove?
params.save_ica = 1; % Save plots and components
params.ica_cor = 0.8; % Cutoff for correlation with EOG/ECG
params.ica_coh = 0.95; % Cutoff for coherence with EOG/ECG

params.corr_threshold = 0.7; % Correlation threshold for badchannel neighbors
params.z_threshold = 20; % Zmax threshold for badchannel and trial detection
params.opm_std_threshold = 5e-12; % Stddev threshold for badtrial detection
params.squid_std_threshold = 2.5e-12; % Stddev threshold for badtrial detection

params.hpi_freq = 33; % HPI coil frequency
params.hpi_gof = 0.7; % Minimum goodness-of-fit for including coil in hpi analysis. JL: CHANGED FROM 0.9

params.trigger_codes = {2};%{1 [3 11] [5 13]}; % combined oddball-nogo and oddball-go
params.trigger_labels = {'Stim_1_cue_onset'};%{'std' 'oddNoGo' 'oddGo'};

params.src_density = '32'; % Sourcemodel density ('4', '8' or '32') = approximate number of sourmes per hemisphere
params.source_fixedori = true; % use fixed orientation sources (along vertex normals); if false: use three orthogonal sources per location
params.noise_cov = 'empty_room'; % noise cov to use; default= ' ' for prestim, alt: 'resting_state', 'empty_room'
params.inv_method = 'mne';
params.numdipoles = 1;
params.peaks = {};
params.peaks{1}.label = 'poststim';
params.peaks{1}.peak_latency = [-0.5 9];

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = {'*by','*bz'};
params.paradigm = 'WMdistractors';

%% Subjects
[subjects, sessions] = getSubjectsAndSessions(base_data_path, false);

if server
    subs_to_run = [find(cellfun(@(x) strcmp(x,'1196'), subjects)) find(cellfun(@(x) strcmp(x,'1206'), subjects)) find(cellfun(@(x) strcmp(x,'1211'), subjects))];
else
    subs_to_run = 1:length(subjects);
end

%% Loop over subjects
for i_sub = subs_to_run
    params.sub = ['sub-' num2str(i_sub,'%02d')];

    %% Paths
    raw_path = fullfile(base_data_path, subjects{i_sub}, sessions{i_sub,1});
    save_path = fullfile(base_save_path, params.sub);
    save_path_mri = fullfile(save_path, 'MRI');
    hpi_path = fullfile(raw_path, 'hedscan');

    if ~exist(save_path, 'dir');     mkdir(save_path);                      end
    if ~exist(fullfile(save_path,'figs'), 'dir'); mkdir(fullfile(save_path,'figs')); end
    if ~exist(save_path_mri, 'dir'); mkdir(save_path_mri);                  end

    %% Files
    opm_file = char(STATIC.subjectData.OPM_FilePath(i_sub));

    tmp = dir(fullfile(raw_path, 'triux/headshape*'));
    aux_file = fullfile(tmp.folder, tmp.name);

    tmp = dir(fullfile(raw_path, 'hedscan', ['*' params.paradigm 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']));
    if ~isempty(tmp)
        squid_file = fullfile(tmp.folder, tmp.name);
    else
        squid_file = [];
    end

    %% Read and preproc - OPM
    params.modality = 'opm';
    params.layout = 'fieldlinebeta2bz_helmet.mat';
    params.chs = {'*by','*bz'};

    if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica.mat']), 'file')
        ft_hastoolbox('mne', 1);

        disp('Reading OPM file...')
        data_epo = read_osMEG(opm_file, {}, save_path, params, ...
            trialinfo(trialinfo.subject==STATIC.subjectID(i_sub),:));

        disp('Running ICA ...')
        if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1
            params.manual_ica = 1;
            params.save_ica = 1;
        end
        data_ica = ica_MEG(data_epo, save_path, params);
        save(fullfile(save_path, [params.paradigm '_data_ica']), 'data_ica', '-v7.3'); disp('done');
        clear data_epo
    else
        data_ica = load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
    end

    if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']), 'file')
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = {'*by','*bz'};
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';

        timelocked = timelock(data_ica, save_path, params);
        save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3');
        clear timelocked
    end
    clear data_ica

    if squid
        %% Read and preproc - SQUID-MAG
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'meg';

        if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat']), 'file')
            ft_hastoolbox('mne', 1);

            disp('Reading SQUID file...')
            data_epo = read_cvMEG(squid_file, params);

            disp('Running ICA ...')
            if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1
                params.manual_ica = 1;
                params.save_ica = 1;
            end
            data_ica = ica_MEG(data_epo, save_path, params);
            save(fullfile(save_path, [params.paradigm '_data_ica_squidmag']), 'data_ica', '-v7.3'); disp('done');
            clear data_epo
        else
            data_ica = load(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat'])).data_ica;
        end

        if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked_squidmag.mat']), 'file')
            params.modality = 'squidmag';
            params.layout = 'neuromag306mag.lay';
            params.chs = 'megmag';
            params.amp_scaler = 1e15;
            params.amp_label = 'B [fT]';
            timelocked = timelock(data_ica, save_path, params);
            save(fullfile(save_path, [params.paradigm '_timelocked_squidmag']), 'timelocked', '-v7.3');
            clear timelocked
        end
        clear data_ica
    end

    %% HPI localization
    ft_hastoolbox('mne', 1);

    if exist(fullfile(save_path, 'opm_trans.mat'), 'file') && overwrite.coreg==false
        disp(['Not overwriting OPM transform for ' params.sub]);
    else
        params.include_chs = {'*_bz'  '*_by'};
        opm_trans = fit_hpi(hpi_path, aux_file, save_path, params);

        opm_timelockedT = load(fullfile(save_path, [params.paradigm '_timelocked.mat'])).timelocked;
        for i = 1:length(params.trigger_labels)
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
            opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
            opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
        end

        % Plot sensor layout vs head/source models
        clear headmodels sourcemodel
        headmodels = load(fullfile(save_path_mri, 'headmodels.mat')).headmodels;
        sourcemodel = load(fullfile(save_path_mri, 'sourcemodel.mat')).sourcemodel;

        h = figure;
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red', 'edgeopacity', 0.5, 'unit', 'cm');
        hold on;
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(opm_timelockedT{1}.grad, 'unit', 'cm')
        hold off;
        title('OPM-MEG')
        view([-140 10])
        saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))
        close all

        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        clear sourcemodel headmodels opm_trans
    end

    %% Dipole fits
    ft_hastoolbox('mne', 1);
    if overwrite.dipole==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
    elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']), 'file')
        headmodel = load(fullfile(save_path_mri, 'headmodels.mat')).headmodels.headmodel_meg;
        mri_resliced = load(fullfile(save_path_mri, 'mri_resliced.mat')).mri_resliced;
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;

        for i_peak = 1:length(params.peaks)
            peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{i_peak}.label])).peak;
            fit_dipoles(save_path, opm_timelockedT, headmodel, mri_resliced, params);
            clear peak_opm
        end
        clear opm_timelockedT
    end

    %% Compute empty room covariance
    er_cov_file = fullfile(save_path, [params.sub '_ER_opm.mat']);
    if ~exist(er_cov_file, 'file')
        er_file = char(STATIC.subjectData.Empty_room_FilePaths(i_sub));
        cfg = [];
        cfg.dataset = er_file;
        cfg.coordsys = 'dewar';
        cfg.coilaccuracy = 0;
        er_raw = ft_preprocessing(cfg);
        cfg = [];
        cfg.channel = {'*bz'};
        cfg.covariance = 'yes';
        er_tl = ft_timelockanalysis(cfg, er_raw);
        opm_ER_cov = er_tl.cov;
        save(er_cov_file, 'opm_ER_cov');
        clear er_raw er_tl opm_ER_cov
    end

    %% MNE
    ft_hastoolbox('mne', 1);
    if exist(fullfile(save_path, 'opm_mne_peaks.mat'), 'file') && overwrite.mne==false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']), 'file')
        clear headmodel sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(save_path_mri, 'sourcemodel.mat')).sourcemodel;
        sourcemodel_inflated = load(fullfile(save_path_mri, 'sourcemodel_inflated.mat')).sourcemodel_inflated;
        headmodel = load(fullfile(save_path_mri, 'headmodels.mat')).headmodels.headmodel_meg;
        sourcemodel.unit = 'cm';
        sourcemodel_inflated.unit = 'cm';

        % OPM
        clear opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;

        for i = 1:length(opm_timelockedT)
            if exist(er_cov_file, 'file')
                opm_timelockedT{i}.cov_ER = load(er_cov_file).opm_ER_cov;
            end
        end

        params.modality = 'opm';
        params.chs = '*bz';
        params.save_mne = true;
        fit_mne(save_path, opm_timelockedT, headmodel, sourcemodel, sourcemodel_inflated, params);

        if squid
            % SQUID
            clear squid_timelocked
            squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;

            for i = 1:length(squid_timelocked)
                squid_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squid_RS_cov;
                if exist(fullfile(save_path, [params.sub '_ER_squid.mat']), 'file')
                    squid_timelocked{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_squid.mat'])).squid_ER_cov;
                end
            end

            params.modality = 'squidgrad';
            params.chs = 'meggrad';
            fit_mne(save_path, squid_timelocked, headmodel, sourcemodel, sourcemodel_inflated, params);
        end
    end
end

% save(fullfile(base_save_path, 'group_results.mat'), 'grp_tag_opm','grp_tag_squid','grp_SNR_opm','grp_SNR_squid','grp_pp_opm','grp_pp_squid');

%% clear and close all, then exit to free memory
close all
clear all
exit

%% Functions
function [subjects, sessions] = getSubjectsAndSessions(folderPath, natmeg)
    if natmeg
        subjectFolders = dir(fullfile(folderPath, 'NatMEG_*'));
    else
        subjectFolders = dir(fullfile(folderPath));
    end
    subjects = {};
    sessions = {};

    i_sub = 0;
    for i = 1:length(subjectFolders)
        if subjectFolders(i).isdir && length(subjectFolders(i).name)>2
            i_sub = i_sub + 1;
            subjects{i_sub,1} = subjectFolders(i).name;

            sessionFolders = dir(fullfile(folderPath, subjectFolders(i).name));
            i_ses = 0;
            for j = 1:length(sessionFolders)
                if sessionFolders(j).isdir && length(sessionFolders(j).name) == 6 && all(isstrprop(sessionFolders(j).name, 'digit'))
                    i_ses = i_ses + 1;
                    sessions{i_sub, i_ses} = sessionFolders(j).name;
                end
            end
        end
    end
end
