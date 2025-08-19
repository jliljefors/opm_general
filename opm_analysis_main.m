%% Reset all and initialize
clear all
close all
restoredefaultpath

STATIC = StaticVarClass_WM_shift_OPM();
addpath(fullfile(STATIC.toolboxpath,"opm_general/"))
% Base paths

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

% Johan
base_data_path = 'G:\SV10_OPM_Distractor_Processing\SV10_RawData';
base_save_path = 'G:\SV10_OPM_Distractor_Processing\SV10_ProcessedData\source_space';
base_matlab_path = STATIC.toolboxpath;
project_scripts_path = 'C:\Github\WM_Shifting_distractors\Matlab\';

% Set up fieldtrip
% addpath(fullfile(base_matlab_path,'fieldtrip')) % Fieldtrip pathr
% addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
% addpath(project_scripts_path)
addpath(fullfile(base_matlab_path,'fieldtrip-20250507/')) % Fieldtrip path
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

% Overwrite
overwrite = [];
overwrite.preproc = true;
overwrite.timelock = true;
overwrite.mri = true;
overwrite.coreg = true;
overwrite.sourcerec = true;

% Params
params = [];
params.pre = 1; % Trial prestim in seconds
params.post = 8.5; % Trial poststim in seconds
params.pad = 0.2; % Trial (pre and post) padding in seconds
params.delay = 0.00; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

params.filter = [];
params.filter.hp_freq = 1; % Highpass cutoff frequency
params.filter.lp_freq = 150; % Lowpass cutoff frequency
params.filter.notch = sort([50 60 100 120]); % Notch (bandstop) filter frequencies

params.apply_hfc = true; % Apply Homogenous Field Correction
params.hfc_order = 2; % Order for Homogenous Field Correction: 1..3

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

params.trigger_codes = 2; % Trigger values to timelock
params.trigger_labels = {'Stim_1'}; % Labels corresponding to the trigger values

params.src_density = '8'; % Sourcemodel density ('4', '8' or '32') = approximate number of sources per hemisphere

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';

% Subjects + dates
%subjects = {'NatMEG_0953'}; % List of subjects to loop through (semicolon separated)
%sessions = {'241104'}; % List of sessions; if multiple per subject define as: {'sub1_ses1' 'sub1_ses2'; 'sub2_ses1' 'sub2_ses2'};
% [subjects, sessions] = getSubjectsAndSessions(base_data_path);

subjects = {'NatMEG_1118';'NatMEG_1212';'NatMEG_1213'};
sessions = {'241122'; '250320';'250325'};

paradigms = {'WM_Shift'}; % Paradigms to analyze for all participants and sessions

if server
    subs_to_run = [find(cellfun(@(x) strcmp(x,'1196'), subjects)) find(cellfun(@(x) strcmp(x,'1206'), subjects)) find(cellfun(@(x) strcmp(x,'1211'), subjects))];
else
    %subs_to_run = find(cellfun(@str2num, subjects)>1082)';
    subs_to_run = find(cellfun(@(x) strcmp(x,'1118'), subjects));
end
ses_cnt = 0;

%% Loop over subjects
for i_sub = 2
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
        % raw_path = fullfile(base_data_path, ['NatMEG_' subjects{i_sub}], sessions{i_sub,i_ses});
        % raw_path = fileparts(STATIC.subjectData.MEGFilePath(i_sub));
        save_path = fullfile(base_save_path, params.sub, params.ses);
        if ~exist(save_path, 'dir')
           mkdir(save_path)
        end
        if ~exist(fullfile(save_path,'figs'), 'dir')
           mkdir(fullfile(save_path,'figs'))
        end
        for i_paradigm = 1:length(paradigms)
            % if server % on server
            %     tmp = dir(fullfile(raw_path,'opm',['*' paradigms{i_paradigm} 'OPM_raw.fif']));
            %     opm_files{i_paradigm} = fullfile(tmp.folder,tmp.name); % opm files 
            %     squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_tsss_mc.fif']); % corresponding aux files containing EOG/ECG
            % else % on laptop
            %     opm_files{i_paradigm} = fullfile(raw_path,'osmeg',[paradigms{i_paradigm} 'OPM_raw.fif']); % opm files 
            %     squid_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']); % corresponding aux files containing EOG/ECG
            % end
            opm_files{i_paradigm} = STATIC.subjectData.OPM_FilePath(i_sub);
            aux_files{i_paradigm} = STATIC.subjectData.MEGFilePath(i_sub);
        end
        hpi_path = STATIC.subjectData.HPI_FilePath(i_sub);
        mri_path = fileparts(STATIC.subjectData.Dicom_Path(i_sub));
        mri_file = STATIC.subjectData.Dicom_Path(i_sub);
        
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
                    cfg = [];
                    cfg.channel = {'EOG*', 'ECG*', 'EEG*'};
                    data_ExG = ft_selectdata(cfg,data_epo);
    
                    cfg = [];
                    cfg.channel         = '*bz';
  	                cfg.order           = params.hfc_order;
                    cfg.updatesens      = 'yes';
                    cfg.residualcheck   = 'no';
                    data_epo = ft_denoise_hfc(cfg,data_epo);
    
                    cfg = [];
                    data_epo = ft_appenddata(cfg,data_epo,data_ExG);
                    clear data_ExG
                end
    
                % AMM
                if params.apply_amm
                    cfg = [];
                    cfg.channel = {'EOG*', 'ECG*', 'EEG*'};
                    data_ExG = ft_selectdata(cfg,data_epo);
    
                    cfg = [];
                    cfg.channel         = '*bz';
                    cfg.updatesens      = 'yes';
                    data_epo = ft_denoise_amm(cfg,data_epo);
    
                    cfg = [];
                    data_epo = ft_appenddata(cfg,data_epo,data_ExG);
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
    
                % %% CAPSI
                % params.modality = 'opmLP';
                % cfg = [];
                % cfg.lpfilter        = 'yes';         
                % cfg.lpfreq          = 20;
                % cfg.demean          = 'yes';
                % data_lp = ft_preprocessing(cfg, data_ica);
                % cfg = [];
                % cfg.latency = [-0.1 0.5];
                % data_lp = ft_selectdata(cfg, data_lp);
                % timelocked_opm_lp = timelock(data_lp, save_path, params); 
                % clear data_lp
                % 
                % params.modality = 'opmHP';
                % cfg = [];
                % cfg.hpfilter        = 'yes';         
                % cfg.hpfreq          = 30;
                % cfg.hpinstabilityfix = 'reduce';
                % cfg.demean          = 'yes';
                % data_hp = ft_preprocessing(cfg, data_ica);
                % timelocked_opm_hp = timelock(data_hp, save_path, params); 
                % clear data_hp
                % 
                % for i_trigger = 1:length(params.trigger_codes)
                %     cfg = [];
                %     cfg.channel = '*bz';
                %     cfg.output = 'pow';
                %     cfg.method = 'mtmfft';
                %     cfg.taper = 'hanning';
                %     cfg.foi = 30:1:50;
                %     cfg.pad = 2;
                %     freq = ft_freqanalysis(cfg, timelocked_opm_hp{i_trigger});   
                %     h = figure;
                %     plot(freq.freq,freq.powspctrm)
                %     xlabel('Frequency (Hz)')
                %     ylabel('Power (T^2)')
                %     [peak_pow, peak_ch] = max(max(freq.powspctrm, [], 2));
                %     noise_pow = mean(freq.powspctrm(peak_ch,[find(freq.freq==32) find(freq.freq==48)]));
                %     snr = peak_pow/noise_pow;
                %     title(['OPM spectrum (peak SNR: ' num2str(snr) ')'])
                %     saveas(h, fullfile(save_path, 'figs', [params.paradigm '_freqTag_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
                %     close all
                %     clear freq
                % end
            end
            clear data_ica
    
            % %% Read and preproc - SQUID-MAG
            % params.modality = 'squidmag';
            % params.layout = 'neuromag306mag.lay';
            % params.chs = 'megmag';
            % 
            % if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat']),'file')
            %     ft_hastoolbox('mne', 1);
            % 
            %     % Read data 
            %     disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(paradigms)) '...'])
            %     data_epo = read_cvMEG(squid_files{i_paradigm}, params); % Read data
            % 
            %     % Reject jump trials
            %     cfg = [];
            %     cfg.channel = {'megmag'};
            %     cfg.metric = 'maxzvalue';
            %     cfg.preproc.medianfilter  = 'yes';
            %     cfg.preproc.medianfiltord  = 9;
            %     cfg.preproc.absdiff       = 'yes';
            %     cfg.threshold = params.z_threshold;
            %     [cfg,badtrl_jump] = ft_badsegment(cfg, data_epo);
            %     data_epo = ft_rejectartifact(cfg,data_epo);
            % 
            %     % Reject noisy trials
            %     cfg = [];
            %     cfg.channel = {'megmag'};
            %     cfg.metric = 'std';
            %     cfg.threshold = params.squid_std_threshold;
            %     [cfg,badtrl_std] = ft_badsegment(cfg, data_epo);
            %     data_epo = ft_rejectartifact(cfg,data_epo);
            % 
            %     % Remove bad trials
            %     [~,idx]=ismember(data_epo.sampleinfo,badtrl_jump,'rows');
            %     badtrl_jump = find(idx);
            %     [~,idx]=ismember(data_epo.sampleinfo,badtrl_std,'rows');
            %     badtrl_std = find(idx);
            %     save(fullfile(save_path, [params.paradigm '_badtrls']), ...
            %         'badtrl_jump', ...
            %         'badtrl_std', "-v7.3"); 
            % 
            %     % ICA
            %     disp('Running ICA ...')
            %     if sum(contains(data_epo.label,'EOG'))<1 || sum(contains(data_epo.label,'ECG'))<1 % No ExG data
            %         params.manual_ica = 1;
            %         params.save_ica = 1;
            %     end
            %     data_ica = ica_MEG(data_epo, save_path, params);
            %     save(fullfile(save_path, [params.paradigm '_data_ica_squidmag']), 'data_ica',"-v7.3"); disp('done');
            %     clear data_epo
            % else
            %     data_ica = load(fullfile(save_path, [params.paradigm '_data_ica_squidmag.mat'])).data_ica;
            % end
            % 
            % if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
            %     params.modality = 'squidmag';
            %     params.layout = 'neuromag306mag.lay';
            %     params.chs = 'megmag';
            %     params.amp_scaler = 1e15;
            %     params.amp_label = 'B [fT]';
            %     timelocked = timelock(data_ica, save_path, params);
            %     save(fullfile(save_path, [params.paradigm '_timelocked_squidmag']), 'timelocked', '-v7.3'); 
            %     clear timelocked
            % 
            %     %% CAPSI
            %     params.modality = 'squidLP';
            %     cfg = [];
            %     cfg.lpfilter        = 'yes';         
            %     cfg.lpfreq          = 20;
            %     cfg.demean          = 'yes';
            %     data_lp = ft_preprocessing(cfg, data_ica);
            %     cfg = [];
            %     cfg.latency = [-0.1 0.5];
            %     data_lp = ft_selectdata(cfg, data_lp);
            %     timelocked_squid_lp = timelock(data_lp, save_path, params); 
            %     clear data_lp
            % 
            %     params.modality = 'squidHP';
            %     cfg = [];
            %     cfg.hpfilter        = 'yes';         
            %     cfg.hpfreq          = 30;
            %     cfg.hpinstabilityfix = 'reduce';
            %     cfg.demean          = 'yes';
            %     data_hp = ft_preprocessing(cfg, data_ica);
            %     timelocked_squid_hp = timelock(data_hp, save_path, params); 
            %     clear data_hp
            % 
            %     for i_trigger = 1:length(params.trigger_codes)
            %         cfg = [];
            %         cfg.channel = 'megmag';
            %         cfg.output = 'pow';
            %         cfg.method = 'mtmfft';
            %         cfg.taper = 'hanning';
            %         cfg.foi = 30:1:50;
            %         cfg.pad = 2;
            %         freq = ft_freqanalysis(cfg, timelocked_squid_hp{i_trigger});          
            %         h = figure;
            %         plot(freq.freq,freq.powspctrm)
            %         xlabel('Frequency (Hz)')
            %         ylabel('Power (T^2)')
            %         [peak_pow, peak_ch] = max(max(freq.powspctrm, [], 2));
            %         noise_pow = mean(freq.powspctrm(peak_ch,[find(freq.freq==32) find(freq.freq==48)]));
            %         snr = peak_pow/noise_pow;
            %         title(['SQUID spectrum (peak SNR: ' num2str(snr) ')'])
            %         saveas(h, fullfile(save_path, 'figs', [params.paradigm '_freqTag_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
            %         close all
            %         clear freq_pre freq_tag
            %     end
            % end
            % clear data_ica
            
        %     %% CAPSI
        %     for i_trigger = 1:length(params.trigger_codes)
        %         % Butterfly comparison
        %         pp_opm = 1e15*max(max(timelocked_opm_lp{i_trigger}.avg))-min(min(timelocked_opm_lp{i_trigger}.avg));
        %         pp_squid = 1e15*max(max(timelocked_squid_lp{i_trigger}.avg))-min(min(timelocked_squid_lp{i_trigger}.avg));
        %         h = figure;
        %         subplot(2,1,1)
        %         plot(timelocked_opm_lp{i_trigger}.time*1e3,timelocked_opm_lp{i_trigger}.avg*1e15)
        %         ylabel('fT')
        %         xlabel('ms')
        %         title(['OPM: ' params.trigger_labels{i_trigger} ' (pk-pk=' num2str(pp_opm) 'fT)'])
        %         ylimits = ylim;
        %         subplot(2,1,2)
        %         plot(timelocked_squid_lp{i_trigger}.time*1e3,timelocked_squid_lp{i_trigger}.avg*1e15)
        %         ylabel('fT')
        %         xlabel('ms')
        %         ylim(ylimits)
        %         title(['SQUID: ' params.trigger_labels{i_trigger} ' (pk-pk=' num2str(pp_squid) 'fT)'])
        %         saveas(h, fullfile(save_path, 'figs', [params.paradigm '_ButterflyComp_trig-' params.trigger_labels{i_trigger} '.jpg']))
        %         close all
        % 
        %         % FFT comparison
        %         cfg = [];
        %         cfg.channel = '*bz';
        %         cfg.output = 'pow';
        %         cfg.method = 'mtmfft';
        %         cfg.taper = 'hanning';
        %         cfg.foi = 30:1:50;
        %         cfg.pad = 2;
        %         freq1 = ft_freqanalysis(cfg, timelocked_opm_hp{i_trigger});
        %         [peak_pow_opm, peak_ch] = max(max(freq1.powspctrm, [], 2));
        %         noise_pow = mean(freq1.powspctrm(peak_ch,[find(freq2.freq<=34) find(freq2.freq>=46)]));
        %         snr_opm = peak_pow_opm/noise_pow;
        %         cfg = [];
        %         cfg.channel = 'megmag';
        %         cfg.output = 'pow';
        %         cfg.method = 'mtmfft';
        %         cfg.taper = 'hanning';
        %         cfg.foi = 30:1:50;
        %         cfg.pad = 2;
        %         freq2 = ft_freqanalysis(cfg, timelocked_squid_hp{i_trigger});  
        %         [peak_pow_squid, peak_ch] = max(max(freq2.powspctrm, [], 2));
        %         noise_pow = mean(freq2.powspctrm(peak_ch,[find(freq2.freq<=34) find(freq2.freq>=46)]));
        %         snr_squid = peak_pow_squid/noise_pow;
        %         h = figure;
        %         subplot(2,1,1)
        %         plot(freq1.freq,freq1.powspctrm)
        %         ylabel('T^2')
        %         xlabel('Hz')
        %         title(['OPM: ' params.trigger_labels{i_trigger} ' (' num2str(snr_opm) ')'])
        %         ylimits = ylim;
        %         subplot(2,1,2)
        %         plot(freq2.freq,freq2.powspctrm)
        %         ylabel('T^2')
        %         xlabel('Hz')
        %         ylim(ylimits)
        %         title(['SQUID: ' params.trigger_labels{i_trigger} ' (' num2str(snr_squid) ')'])
        %         saveas(h, fullfile(save_path, 'figs', [params.paradigm '_FreqTagComp_trig-' params.trigger_labels{i_trigger} '.jpg']))
        %         close all
        %         clear freq1 freq2
        % 
        %         grp_pp_squid(ses_cnt,i_trigger) = pp_squid;
        %         grp_pp_opm(ses_cnt,i_trigger) = pp_opm;
        %         grp_SNR_squid(ses_cnt,i_trigger) = snr_squid;
        %         grp_SNR_opm(ses_cnt,i_trigger) = snr_opm;
        %         grp_tag_squid(ses_cnt,i_trigger) = peak_pow_squid;
        %         grp_tag_opm(ses_cnt,i_trigger) = peak_pow_opm;
        %     end
        %     clear timelocked_opm_lp timelocked_opm_hp timelocked_squid_lp timelocked_squid_hp
        % end



        %% MRI
        % Assumes you have a brain sourcemodel (https://www.fieldtriptoolbox.org/workshop/practicalmeeg2022/handson_anatomy/) and aligns that to headspace as well. 
        
        % Plot standard freesurfer output
        freesurfer_output_path = fullfile(STATIC.wsl_path,string(STATIC.subjectID(i_sub)),string(STATIC.subjectID(i_sub)),"surf");
      
        ft_hastoolbox('mne',1);
        if overwrite.mri==true
            ft_hastoolbox('mne', 1);
            prepare_mri(mri_file,opm_files{1},aux_files{1},save_path)   % REMOVE PLOTTING IN prepare_mri as thats what uses .grad
        end

        %%  FREESURFER
        % Note this script is best run on the raw dicoms eg can be done as soon as the MRI is completed and is not dependend on fieldtrip
        if ~exist( fullfile(STATIC.wsl_path,string(STATIC.subjectID(i_sub))),'dir')
            mkdir(fullfile(STATIC.wsl_path,string(STATIC.subjectID(i_sub))));
        end
        copyfile(fullfile(fileparts(save_path),'freesurfer\sub02.mgh'),... 
            fullfile(STATIC.wsl_path,string(STATIC.subjectID(i_sub))))
        % Copy nifti file to WSL directory, start WSL and run script

        %% Sourcemodel
            % Read and transform cortical restrained source model
            wsl_path = fullfile(STATIC.wsl_path,string(STATIC.subjectID(i_sub)),string(STATIC.subjectID(i_sub)));
            files = dir(fullfile(wsl_path,'workbench'));
            for i = 1:length(files)
                if endsWith(files(i).name,['.L.midthickness.' params.src_density 'k_fs_LR.surf.gii'])
                    fn = fullfile(wsl_path,'workbench',files(i).name);
                    [~,filename,ext] = fileparts(fn);
                    filename = fullfile(getenv('USERPROFILE'),'Downloads',strcat(filename,ext));
                    copyfile(fn,filename)
                    copyfile(strrep(fn,'.L.','.R.'),strrep(filename,'.L.','.R.'))
                end
                if endsWith(files(i).name,['.L.aparc.' params.src_density 'k_fs_LR.label.gii'])
                    fn = fullfile(wsl_path,'workbench',files(i).name);
                    [~,filename2,ext] = fileparts(fn);
                    filename2 = fullfile(getenv('USERPROFILE'),'Downloads',strcat(filename2,ext));
                    copyfile(fn,filename2)
                    copyfile(strrep(fn,'.L.','.R.'),strrep(filename2,'.L.','.R.'))
                end
            end
            clear sourcemodel mri_resliced
            mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
            sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

            aparc_L = ft_read_atlas({num2str(filename2),num2str(filename)});
            aparc_R = ft_read_atlas({num2str(strrep(filename2,'.L.','.R.')),num2str(strrep(filename,'.L.','.R.'))});
            tmp = ft_read_atlas(num2str(strrep(filename2, '.L.', '.R.')),'format','caret_label');
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
        
            % SINCE THE MRI I HAVE PASSED TO FREESURFER IS ALREADY ALIGNED, I SHOULDNT HAVE TO DO THIS
            % BUT NEED TO PLOT RESULT!
            % T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
            % sourcemodel = ft_transform_geometry(T, sourcemodel);
            % sourcemodel.inside = true(size(sourcemodel.pos,1),1);
            save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
    
            clear mri_resliced sourcemodel T atlas tmp aparc_L aparc_R filename filename2
        end
        
        %% HPI localization
        ft_hastoolbox('mne',1);
        
        if overwrite.coreg==true 
            ft_hastoolbox('mne', 1);
            % params.include_chs = load(fullfile(save_path, ['include_chs' num2str(length(opm_files))])).include_chs;
            data_ica_ds=load(fullfile(save_path, 'WM_Shift_data_ica.mat')).data_ica;
            params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));

            clear data_ica
            opm_trans = fit_hpi(hpi_path, aux_files{1}, save_path, params);
            % 
            % for i_file = 1:length(opm_files)
            %     clear timelocked
            %     data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_auditory' num2str(i_file) '.mat'])).data_ica;
            %     data_ica.grad.chanpos = opm_trans.transformPointsForward(data_ica.grad.chanpos);
            %     data_ica.grad.coilpos = opm_trans.transformPointsForward(data_ica.grad.coilpos);
            %     data_ica.grad.chanori = (opm_trans.Rotation'*data_ica.grad.chanori')';
            %     data_ica.grad.coilori = (opm_trans.Rotation'*data_ica.grad.coilori')';
            %     save(fullfile(save_path, [params.sub '_' params.modality '_auditoryT' num2str(i_file)]),'data_ica',"-v7.3");
            %     cfg = [];
            %     cfg.covariance = 'yes';
            %     cfg.covariancewindow = [-params.pre 0];
            %     timelocked = ft_timelockanalysis(cfg,data_ica);
            %     save(fullfile(save_path, [params.sub '_' params.modality '_auditory_timelockedT' num2str(i_file)]),'timelocked',"-v7.3");
            %     clear data_ica timelocked
            % end
    
            % Plot source and head models
            clear headmodels sourcemodel
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;
        
            timelocked = load(fullfile(save_path,"timelocked.mat")).timelocked;
           
            timelocked.grad = ft_transform_geometry(opm_trans.A, timelocked.grad);    % JL: NEED TO APPLY THE TRANSFORM ABOVE TO ALIGN THE DATA. 
            save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3');
            
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
            sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            pos = sourcemodel.pos;
            tri = sourcemodel.tri;
            roi = sourcemodel.brainstructurelabel;
            
            N_rois = length(sourcemodel.brainstructurelabel);
            N_sources = size(pos, 1);
            mapping_matrix = zeros(N_rois, N_sources);
            for i = 1:N_rois
                mapping_matrix(i,sourcemodel.brainstructure==i) = 1;
            end
            roi_counts = sum(mapping_matrix, 2);
            mapping_matrix = mapping_matrix ./ repmat(roi_counts,[1 size(mapping_matrix,2)]);
        
            for i_file = 1:length(opm_files)
                if i_ses == 1
                    data_set = '';
                else
                    data_set = num2str(i_file);
                end
            
                % timelocked = load(fullfile(save_path, [params.sub '_' params.modality '_motorimag_timelockedT' data_set '.mat'])).timelocked;
                timelocked = load(fullfile(save_path,"timelocked.mat")).timelocked; % JL: LETS TRY THIS
                % mne = fit_mne_opmbci(timelocked,headmodels,sourcemodel,params); 
                mne = fit_mne(timelocked,headmodels,sourcemodel,params); 
                mne_inv = zeros(length(mne.avg.filter),length(mne.avg.label));
                for i = 1:length(mne.avg.filter)
                    mne_inv(i,:) = mne.avg.filter{i};
                end
                % Plot timelocked
                h = figure ;
                plot(timelocked.time,timelocked.avg)
                title('Motor imagery timelocked')
                saveas(h, fullfile(save_path, 'figs', 'opm_motorimag_timelocked.jpg'))
                close all
                clear mne timelocked
                data_ica = load(fullfile(save_path, [params.sub '_' params.modality '_motorimag' data_set '.mat'])).data_ica;
                trial = cell(length(data_ica.trial),1);
                for i_trl = 1:length(data_ica.trial)
                    trial{i_trl} = mapping_matrix*((mne_inv*data_ica.trial{i_trl}));
                end
                time = data_ica.time{1};
                label = params.trigger_labels(data_ica.trialinfo-64);
                save(fullfile(save_path, ['motorimag_mne' num2str(i_file) '.mat']),'mne_inv','trial','time','label','pos','tri','roi','-v7.3');
                clear mne_inv trial time label  data_ica
            end
            clear pos tri roi sourcemodel headmodels mapping_matrix
        end
        close all
    end
end

% save(fullfile(base_save_path, 'group_results.mat'), 'grp_tag_opm','grp_tag_squid','grp_SNR_opm','grp_SNR_squid','grp_pp_opm','grp_pp_squid');

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
