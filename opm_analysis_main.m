%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    server = true;
    % Server:
    base_data_path = '/archive/23108_CAPSI/MEG/';
    base_save_path = '/home/share/capsi_share';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/opm_general';
else
    server = false;
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/21099_opm/MEG';
    base_save_path = '/Users/christophpfeiffer/data_local/24110_opm_auditory';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_general';
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Overwrite
overwrite = [];
overwrite.preproc = true;
overwrite.timelock = true;
overwrite.mri = false;
overwrite.coreg = false;
overwrite.sourcerec = false;

%% Params
params = [];
params.pre = 0.2; % Trial prestim in seconds
params.post = 0.8; % Trial poststim in seconds
params.pad = 0.2; % Trial (pre and post) padding in seconds
params.delay = 0.01; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

params.filter = [];
params.filter.hp_freq = 1; % Highpass cutoff frequency
params.filter.lp_freq = 50; % Lowpass cutoff frequency
params.filter.notch = sort([50 60 100 120]); % Notch (bandstop) filter frequencies

params.apply_hfc = true; % Apply Homogenous Field Correction
params.hfc_order = 1; % Order for Homogenous Field Correction: 1..3

params.apply_amm = false; % Apply Adaptive Multipole Models

params.n_comp = 40; % Number of ICA components 
params.manual_ica = false; % Manually select ICA components to remove?
params.ica_cor = 0.8; % Cutoff for correlation with EOG/ECG 
params.ica_coh = 0.95; % Cutoff for coherence with EOG/ECG 
params.save_ica = 1; % Save plots and components

params.corr_threshold = 0.7; % Correlation threshold for badchannel neighbors
params.z_threshold = 20; % Zmax threshold for badchannel and trial detection
params.opm_std_threshold = 5e-12; % Stddev threshold for badtrial detection
params.squid_std_threshold = 2.5e-12; % Stddev threshold for badtrial detection

params.hpi_freq = 33; % HPI coil frequency
params.hpi_gof = 0.9; % Minimum goodness-of-fit for including coil in hpi analysis

params.trigger_codes = [1 3 5 11 13]; % Trigger values to timelock
params.trigger_labels = {'std', 'ngl', 'gol', 'ngh', 'goh'}; % Labels corresponding to the trigger values

params.src_density = '8'; % Sourcemodel density ('4', '8' or '32') = approximate number of sources per hemisphere

params.cov = 'resting_state';

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';

%% Subjects + dates
%subjects = {'NatMEG_0953'}; % List of subjects to loop through (semicolon separated)
%sessions = {'241104'}; % List of sessions; if multiple per subject define as: {'sub1_ses1' 'sub1_ses2'; 'sub2_ses1' 'sub2_ses2'};
[subjects, sessions] = getSubjectsAndSessions(base_data_path);

paradigms = {'AudOdd'}; % Paradigms to analyze for all participants and sessions

if server
    subs_to_run = [find(cellfun(@(x) strcmp(x,'1196'), subjects)) find(cellfun(@(x) strcmp(x,'1206'), subjects)) find(cellfun(@(x) strcmp(x,'1211'), subjects))];
else
    %subs_to_run = find(cellfun(@str2num, subjects)>1082)';
    subs_to_run = find(cellfun(@(x) strcmp(x,'0953'), subjects));
end
ses_cnt = 0;

%% Loop over subjects
for i_sub = subs_to_run
    % Loop over sessions
    for i_ses = 1:length(sessions(i_sub,:))
        if isempty(sessions{i_sub,i_ses})
            disp(['No session defined! Skipping sub-' num2str(i_sub,'%02d') '_ses-' num2str(i_ses,'%02d')])
            continue % Skip iteration if no session defined
        end
    
        ses_cnt = ses_cnt + 1;

        %% Loop over subjects
        params.sub = ['sub-' num2str(i_sub,'%02d')];
        params.ses = ['ses-' num2str(i_ses,'%02d')];
        
        %% Paths
        raw_path = fullfile(base_data_path, ['NatMEG_' subjects{i_sub}], sessions{i_sub,i_ses});
        save_path = fullfile(base_save_path, params.sub, params.ses);
        if ~exist(save_path, 'dir')
           mkdir(save_path)
        end
        if ~exist(fullfile(save_path,'figs'), 'dir')
           mkdir(fullfile(save_path,'figs'))
        end
        for i_paradigm = 1:length(paradigms)
            if server % on server
                tmp = dir(fullfile(raw_path,'opm',['*' paradigms{i_paradigm} 'OPM_raw.fif']));
                opm_files{i_paradigm} = fullfile(tmp.folder,tmp.name); % opm files 
                squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_tsss_mc.fif']); % corresponding aux files containing EOG/ECG
            else % on laptop
                opm_files{i_paradigm} = fullfile(raw_path,'osmeg',[paradigms{i_paradigm} 'OPM_raw.fif']); % opm files 
                squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']); % corresponding aux files containing EOG/ECG
            end
            aux_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'EEG.fif']); % corresponding aux files containing EOG/ECG
        end
        hpi_path = fullfile(raw_path,'osmeg');
        mri_path = '/Volumes/dataarchvie/23106_opmbci/NatMEG_0953/mri/';
        mri_file = fullfile(mri_path,'mri','orig','001.mgz');
        
        for i_paradigm = 1:length(paradigms)
            params.paradigm = paradigms{i_paradigm};
            
            %% Read and preproc
            params.modality = 'opm';
            params.layout = 'fieldlinebeta2bz_helmet.mat';
            params.chs = '*bz';
            
            if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica.mat']),'file')
                ft_hastoolbox('mne', 1);
    
                % Read data 
                disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigms)) '...'])
                [data_epo, badchs_opm] = read_osMEG(opm_files{i_paradigm}, aux_files{i_paradigm}, save_path, params); % Read data
                
                % Reject bad channels
                cfg = [];
                cfg.channel = setdiff(data_epo.label,badchs_opm);
                data_epo = ft_selectdata(cfg, data_epo);         
        
                % HFC
                if params.apply_hfc
                    % Save ECG, EOG and EEG channels
                    cfg = [];
                    cfg.channel = {'EOG*', 'ECG*', 'EEG*'};
                    data_ExG = ft_selectdata(cfg,data_epo);
    
                    cfg = [];
                    cfg.channel         = '*bz';
  	                cfg.order           = params.hfc_order;
                    cfg.updatesens      = 'yes';
                    cfg.residualcheck   = 'no';
                    data_epo = ft_denoise_hfc(cfg,data_epo);
    
                    % Replace ECG, EOG and EEG channels
                    data_epo.label = vertcat(data_epo.label,data_ExG.label);
                    for i = 1:length(data_epo.trial)
                        data_epo.trial{i} = vertcat(data_epo.trial{i}, data_ExG.trial{i}); 
                    end
                    clear data_ExG
                end
    
                % AMM
                if params.apply_amm
                    % Save ECG, EOG and EEG channels
                    cfg = [];
                    cfg.channel = {'EOG*', 'ECG*', 'EEG*'};
                    data_ExG = ft_selectdata(cfg,data_epo);
    
                    cfg = [];
                    cfg.channel         = '*bz';
                    cfg.updatesens      = 'yes';
                    data_epo = ft_denoise_amm(cfg,data_epo);
    
                    % Replace ECG, EOG and EEG channels
                    data_epo.label = vertcat(data_epo.label,data_ExG.label);
                    for i = 1:length(data_epo.trial)
                        data_epo.trial{i} = vertcat(data_epo.trial{i}, data_ExG.trial{i}); 
                    end
                    clear data_ExG
                end
    
                % Reject jump trials
                cfg = [];
                cfg.channel = {'*bz'};
                cfg.metric = 'maxzvalue';
                cfg.preproc.medianfilter  = 'yes';
                cfg.preproc.medianfiltord  = 9;
                cfg.preproc.absdiff       = 'yes';
                cfg.threshold = params.z_threshold;
                [cfg,badtrl_jump] = ft_badsegment(cfg, data_epo);
                data_epo = ft_rejectartifact(cfg,data_epo);
                
                % Reject noisy trials
                cfg = [];
                cfg.channel = {'*bz'};
                cfg.metric = 'std';
                cfg.threshold = params.opm_std_threshold;
                [cfg,badtrl_std] = ft_badsegment(cfg, data_epo);
                data_epo = ft_rejectartifact(cfg,data_epo);
    
                % Remove bad trials
                [~,idx]=ismember(data_epo.sampleinfo,badtrl_jump,'rows');
                badtrl_jump = find(idx);
                [~,idx]=ismember(data_epo.sampleinfo,badtrl_std,'rows');
                badtrl_std = find(idx);
                save(fullfile(save_path, [params.paradigm '_badtrls']), ...
                    'badtrl_jump', ...
                    'badtrl_std', "-v7.3"); 
                
                % ICA
                disp('Running ICA ...')
                if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
                    params.manual_ica = 1;
                    params.save_ica = 1;
                end
                data_ica = ica_MEG(data_epo, save_path, params);
                save(fullfile(save_path, [params.paradigm '_data_ica']), 'data_ica',"-v7.3"); disp('done');
                clear data_epo
            else
                data_ica = load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
            end
            
            if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
                params.modality = 'opm';
                params.layout = 'fieldlinebeta2bz_helmet.mat';
                params.chs = '*bz';
                params.amp_scaler = 1e15;
                params.amp_label = 'B [fT]';
                timelocked = timelock(data_ica, save_path, params);
                save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3'); 
                clear timelocked
    
                %% CAPSI
                params.modality = 'opmLP';
                cfg = [];
                cfg.lpfilter        = 'yes';         
                cfg.lpfreq          = 20;
                cfg.demean          = 'yes';
                data_lp = ft_preprocessing(cfg, data_ica);
                cfg = [];
                cfg.latency = [-0.1 0.5];
                data_lp = ft_selectdata(cfg, data_lp);
                timelocked_opm_lp = timelock(data_lp, save_path, params); 
                clear data_lp
    
                params.modality = 'opmHP';
                cfg = [];
                cfg.hpfilter        = 'yes';         
                cfg.hpfreq          = 30;
                cfg.hpinstabilityfix = 'reduce';
                cfg.demean          = 'yes';
                data_hp = ft_preprocessing(cfg, data_ica);
                timelocked_opm_hp = timelock(data_hp, save_path, params); 
                clear data_hp
    
                for i_trigger = 1:length(params.trigger_codes)
                    cfg = [];
                    cfg.channel = '*bz';
                    cfg.output = 'pow';
                    cfg.method = 'mtmfft';
                    cfg.taper = 'hanning';
                    cfg.foi = 30:1:50;
                    cfg.pad = 2;
                    freq = ft_freqanalysis(cfg, timelocked_opm_hp{i_trigger});   
                    h = figure;
                    plot(freq.freq,freq.powspctrm)
                    xlabel('Frequency (Hz)')
                    ylabel('Power (T^2)')
                    [peak_pow, peak_ch] = max(max(freq.powspctrm, [], 2));
                    noise_pow = mean(freq.powspctrm(peak_ch,[find(freq.freq==32) find(freq.freq==48)]));
                    snr = peak_pow/noise_pow;
                    title(['OPM spectrum (peak SNR: ' num2str(snr) ')'])
                    saveas(h, fullfile(save_path, 'figs', [params.paradigm '_freqTag_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
                    close all
                    clear freq
                end
            end
            clear data_ica
    
            %% Read and preproc - SQUID-MAG
            params.modality = 'squidmag';
            params.layout = 'neuromag306mag.lay';
            params.chs = 'megmag';
    
            if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat']),'file')
                ft_hastoolbox('mne', 1);
    
                % Read data 
                disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigms)) '...'])
                data_epo = read_cvMEG(squid_files{i_paradigm}, params); % Read data
    
                % Reject jump trials
                cfg = [];
                cfg.channel = {'megmag'};
                cfg.metric = 'maxzvalue';
                cfg.preproc.medianfilter  = 'yes';
                cfg.preproc.medianfiltord  = 9;
                cfg.preproc.absdiff       = 'yes';
                cfg.threshold = params.z_threshold;
                [cfg,badtrl_jump] = ft_badsegment(cfg, data_epo);
                data_epo = ft_rejectartifact(cfg,data_epo);
                
                % Reject noisy trials
                cfg = [];
                cfg.channel = {'megmag'};
                cfg.metric = 'std';
                cfg.threshold = params.squid_std_threshold;
                [cfg,badtrl_std] = ft_badsegment(cfg, data_epo);
                data_epo = ft_rejectartifact(cfg,data_epo);
    
                % Remove bad trials
                [~,idx]=ismember(data_epo.sampleinfo,badtrl_jump,'rows');
                badtrl_jump = find(idx);
                [~,idx]=ismember(data_epo.sampleinfo,badtrl_std,'rows');
                badtrl_std = find(idx);
                save(fullfile(save_path, [params.paradigm '_badtrls']), ...
                    'badtrl_jump', ...
                    'badtrl_std', "-v7.3"); 
                
                % ICA
                disp('Running ICA ...')
                if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
                    params.manual_ica = 1;
                    params.save_ica = 1;
                end
                data_ica = ica_MEG(data_epo, save_path, params);
                save(fullfile(save_path, [params.paradigm '_data_ica_squidmag']), 'data_ica',"-v7.3"); disp('done');
                clear data_epo
            else
                data_ica = load(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat'])).data_ica;
            end
            
            if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
                params.modality = 'squidmag';
                params.layout = 'neuromag306mag.lay';
                params.chs = 'megmag';
                params.amp_scaler = 1e15;
                params.amp_label = 'B [fT]';
                timelocked = timelock(data_ica, save_path, params);
                save(fullfile(save_path, [params.paradigm '_timelocked_squidmag']), 'timelocked', '-v7.3'); 
                clear timelocked
    
                %% CAPSI
                params.modality = 'squidLP';
                cfg = [];
                cfg.lpfilter        = 'yes';         
                cfg.lpfreq          = 20;
                cfg.demean          = 'yes';
                data_lp = ft_preprocessing(cfg, data_ica);
                cfg = [];
                cfg.latency = [-0.1 0.5];
                data_lp = ft_selectdata(cfg, data_lp);
                timelocked_squid_lp = timelock(data_lp, save_path, params); 
                clear data_lp
    
                params.modality = 'squidHP';
                cfg = [];
                cfg.hpfilter        = 'yes';         
                cfg.hpfreq          = 30;
                cfg.hpinstabilityfix = 'reduce';
                cfg.demean          = 'yes';
                data_hp = ft_preprocessing(cfg, data_ica);
                timelocked_squid_hp = timelock(data_hp, save_path, params); 
                clear data_hp
    
                for i_trigger = 1:length(params.trigger_codes)
                    cfg = [];
                    cfg.channel = 'megmag';
                    cfg.output = 'pow';
                    cfg.method = 'mtmfft';
                    cfg.taper = 'hanning';
                    cfg.foi = 30:1:50;
                    cfg.pad = 2;
                    freq = ft_freqanalysis(cfg, timelocked_squid_hp{i_trigger});          
                    h = figure;
                    plot(freq.freq,freq.powspctrm)
                    xlabel('Frequency (Hz)')
                    ylabel('Power (T^2)')
                    [peak_pow, peak_ch] = max(max(freq.powspctrm, [], 2));
                    noise_pow = mean(freq.powspctrm(peak_ch,[find(freq.freq==32) find(freq.freq==48)]));
                    snr = peak_pow/noise_pow;
                    title(['SQUID spectrum (peak SNR: ' num2str(snr) ')'])
                    saveas(h, fullfile(save_path, 'figs', [params.paradigm '_freqTag_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
                    close all
                    clear freq_pre freq_tag
                end
            end
            clear data_ica
            
            %% CAPSI
            for i_trigger = 1:length(params.trigger_codes)
                % Butterfly comparison
                pp_opm = 1e15*max(max(timelocked_opm_lp{i_trigger}.avg))-min(min(timelocked_opm_lp{i_trigger}.avg));
                pp_squid = 1e15*max(max(timelocked_squid_lp{i_trigger}.avg))-min(min(timelocked_squid_lp{i_trigger}.avg));
                h = figure;
                subplot(2,1,1)
                plot(timelocked_opm_lp{i_trigger}.time*1e3,timelocked_opm_lp{i_trigger}.avg*1e15)
                ylabel('fT')
                xlabel('ms')
                title(['OPM: ' params.trigger_labels{i_trigger} ' (pk-pk=' num2str(pp_opm) 'fT)'])
                ylimits = ylim;
                subplot(2,1,2)
                plot(timelocked_squid_lp{i_trigger}.time*1e3,timelocked_squid_lp{i_trigger}.avg*1e15)
                ylabel('fT')
                xlabel('ms')
                ylim(ylimits)
                title(['SQUID: ' params.trigger_labels{i_trigger} ' (pk-pk=' num2str(pp_squid) 'fT)'])
                saveas(h, fullfile(save_path, 'figs', [params.paradigm '_ButterflyComp_trig-' params.trigger_labels{i_trigger} '.jpg']))
                close all

                % FFT comparison
                cfg = [];
                cfg.channel = '*bz';
                cfg.output = 'pow';
                cfg.method = 'mtmfft';
                cfg.taper = 'hanning';
                cfg.foi = 30:1:50;
                cfg.pad = 2;
                freq1 = ft_freqanalysis(cfg, timelocked_opm_hp{i_trigger});
                [peak_pow_opm, peak_ch] = max(max(freq1.powspctrm, [], 2));
                noise_pow = mean(freq1.powspctrm(peak_ch,[find(freq2.freq<=34) find(freq2.freq>=46)]));
                snr_opm = peak_pow_opm/noise_pow;
                cfg = [];
                cfg.channel = 'megmag';
                cfg.output = 'pow';
                cfg.method = 'mtmfft';
                cfg.taper = 'hanning';
                cfg.foi = 30:1:50;
                cfg.pad = 2;
                freq2 = ft_freqanalysis(cfg, timelocked_squid_hp{i_trigger});  
                [peak_pow_squid, peak_ch] = max(max(freq2.powspctrm, [], 2));
                noise_pow = mean(freq2.powspctrm(peak_ch,[find(freq2.freq<=34) find(freq2.freq>=46)]));
                snr_squid = peak_pow_squid/noise_pow;
                h = figure;
                subplot(2,1,1)
                plot(freq1.freq,freq1.powspctrm)
                ylabel('T^2')
                xlabel('Hz')
                title(['OPM: ' params.trigger_labels{i_trigger} ' (' num2str(snr_opm) ')'])
                ylimits = ylim;
                subplot(2,1,2)
                plot(freq2.freq,freq2.powspctrm)
                ylabel('T^2')
                xlabel('Hz')
                ylim(ylimits)
                title(['SQUID: ' params.trigger_labels{i_trigger} ' (' num2str(snr_squid) ')'])
                saveas(h, fullfile(save_path, 'figs', [params.paradigm '_FreqTagComp_trig-' params.trigger_labels{i_trigger} '.jpg']))
                close all
                clear freq1 freq2

                grp_pp_squid(ses_cnt,i_trigger) = pp_squid;
                grp_pp_opm(ses_cnt,i_trigger) = pp_opm;
                grp_SNR_squid(ses_cnt,i_trigger) = snr_squid;
                grp_SNR_opm(ses_cnt,i_trigger) = snr_opm;
                grp_tag_squid(ses_cnt,i_trigger) = peak_pow_squid;
                grp_tag_opm(ses_cnt,i_trigger) = peak_pow_opm;
            end
            clear timelocked_opm_lp timelocked_opm_hp timelocked_squid_lp timelocked_squid_hp
        end
        
        %% MRI
        ft_hastoolbox('mne',1);
        if overwrite.mri==true
            ft_hastoolbox('mne', 1);
            prepare_mri(mri_file,aux_files(1),save_path)
    
            % Read and transform cortical restrained source model
            files = dir(fullfile(mri_path,'workbench'));
            for i = 1:length(files)
                if endsWith(files(i).name,['.L.midthickness.' params.src_density 'k_fs_LR.surf.gii'])
                    filename = fullfile(mri_path,'workbench',files(i).name);
                end
                if endsWith(files(i).name,['.L.aparc.' params.src_density 'k_fs_LR.label.gii'])
                    filename2 = fullfile(mri_path,'workbench',files(i).name);
                end
            end
            clear sourcemodel mri_resliced
            mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
            sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
            
            aparc_L = ft_read_atlas({filename2,filename});
            aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
            tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
            n_labels = length(aparc_L.parcellationlabel);
            atlas = [];
            atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
            atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
            atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
            n_labels = length(atlas.parcellationlabel);
            atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
            sourcemodel.brainstructure = atlas.parcellation;
            sourcemodel.brainstructurelabel = atlas.parcellationlabel;
            sourcemodel.brainstructurecolor = atlas.rgba;
        
            T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
            sourcemodel = ft_transform_geometry(T, sourcemodel);
            sourcemodel.inside = true(size(sourcemodel.pos,1),1);
            save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
    
            clear mri_resliced sourcemodel T atlas tmp aparc_L aparc_R filename filename2
        end
        
        %% HPI localization
        ft_hastoolbox('mne',1);
        
        if overwrite.coreg==true 
            ft_hastoolbox('mne', 1);
            params.include_chs = load(fullfile(save_path, ['include_chs' num2str(length(opm_files))])).include_chs;
            clear data_ica
            opm_trans = fit_hpi(hpi_path, aux_files{1}, save_path, params);
        
            for i_file = 1:length(opm_files)
                clear timelocked
                data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_auditory' num2str(i_file) '.mat'])).data_ica;
                data_ica.grad.chanpos = opm_trans.transformPointsForward(data_ica.grad.chanpos);
                data_ica.grad.coilpos = opm_trans.transformPointsForward(data_ica.grad.coilpos);
                data_ica.grad.chanori = (opm_trans.Rotation'*data_ica.grad.chanori')';
                data_ica.grad.coilori = (opm_trans.Rotation'*data_ica.grad.coilori')';
                save(fullfile(save_path, [params.sub '_' params.modality '_auditoryT' num2str(i_file)]),'data_ica',"-v7.3");
                cfg = [];
                cfg.covariance = 'yes';
                cfg.covariancewindow = [-params.pre 0];
                timelocked = ft_timelockanalysis(cfg,data_ica);
                save(fullfile(save_path, [params.sub '_' params.modality '_auditory_timelockedT' num2str(i_file)]),'timelocked',"-v7.3");
                clear data_ica timelocked
            end
    
            % Plot source and head models
            clear headmodels sourcemodel
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            sourcemodel = load(fullfile(save_path, 'sourcemodel.mat')).sourcemodel;
        
            h=figure; 
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
            hold on; 
            ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(timelocked.grad,'unit','cm')
            hold off;
            title('OPM-MEG')
            view([-140 10])
            saveas(h, fullfile(save_path, 'figs', 'opm_layout.jpg'))
            close all
        
            clear timelocked sourcemodel headmodels opm_trans
        end
        
        %% MNE
        ft_hastoolbox('mne',1);
        if overwrite.sourcerec==true
            clear headmodels sourcemodel
            sourcemodel = load(fullfile(save_path, 'sourcemodel')).sourcemodel;
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            
            N_rois = length(sourcemodel.brainstructurelabel);
            N_sources = size(sourcemodel.pos, 1);
            mapping_matrix = zeros(N_rois, N_sources);
            for i = 1:N_rois
                mapping_matrix(i,sourcemodel.brainstructure==i) = 1;
            end
            roi_counts = sum(mapping_matrix, 2);
            mapping_matrix = mapping_matrix ./ repmat(roi_counts,[1 size(mapping_matrix,2)]);
    
            opm_timelocked = load(fullfile(save_path, [params.sub '_opm_auditory_timelockedT.mat'])).timelocked;
            squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_auditory_timelocked.mat'])).timelocked;
            squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_auditory_timelocked.mat'])).timelocked;
            mne = fit_mne(save_path,squidmag_timelocked,squidgrad_timelocked,opm_timelocked,headmodels,sourcemodel,[],params); 
            close all
        end
    end
end

save(fullfile(base_save_path, 'group_results.mat'), 'grp_tag_opm','grp_tag_squid','grp_SNR_opm','grp_SNR_squid','grp_pp_opm','grp_pp_squid');

%% clear and close all, then exit to free memory
close all
clear all
exit

%% Functions
function [subjects, sessions] = getSubjectsAndSessions(folderPath)
    % Get a list of all subject folders
    subjectFolders = dir(fullfile(folderPath, 'NatMEG_*'));
    
    % Initialize cell arrays for subjects and sessions
    subjects = {};
    sessions = {};
    
    i_sub = 0;
    for i = 1:length(subjectFolders)
        if subjectFolders(i).isdir
            i_sub = i_sub + 1;
            subjectID = subjectFolders(i).name(end-3:end);
            subjects{i_sub,1} = subjectID;
            
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
