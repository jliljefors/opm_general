%% Reset all and initialize
clear all
close all
restoredefaultpath

STATIC = StaticVarClass_WM_shift_OPM();
addpath(fullfile(STATIC.toolboxpath,"opm_general/"))
% Base paths

% Johan
base_save_path = 'G:\SV10_OPM_Distractor_Processing\SV10_ProcessedData\source_space';
base_matlab_path = STATIC.toolboxpath;

% Set up fieldtrip

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

subs_to_run = find(cellfun(@(x) strcmp(x,'1118'), subjects));

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
        save_path = fullfile(base_save_path, params.sub, params.ses);
        if ~exist(save_path, 'dir')
            mkdir(save_path)
        end
        if ~exist(fullfile(save_path,'figs'), 'dir')
            mkdir(fullfile(save_path,'figs'))
        end
        for i_paradigm = 1:length(paradigms)
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
                grad = data_epo.grad;
                % Reject bad channels
                cfg = [];
                cfg.channel = setdiff(data_epo.label,badchs_opm);
                params.include_chs = cfg.channel;
                data_epo = ft_selectdata(cfg, data_epo);

                % HFC JL: For some reason the grad structure disappears here
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
                data_ica.grad = grad;  % JL: % Not sure why the .grad structure isnt there, it disappears in the HFC

                save(fullfile(save_path, [params.paradigm '_data_ica']), 'data_ica',"-v7.3"); disp('done');
                clear data_epo
            else
                data_ica = load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
            end

            if overwrite.timelock == true || ~exist(fullfile(save_path, [params.paradigm '_timelocked.mat']),'file')
                params.modality = 'opm';
                params.layout = 'fieldlinebeta2bz_helmet.mat';
                params.chs = '*bz';
                params.chs = [{'*bz'} strcat('-', badchs_opm')];  % JL EXCLUDE BAD CHANNELS
                params.amp_scaler = 1e15;
                params.amp_label = 'B [fT]';
                timelocked = timelock(data_ica, save_path, params); timelocked = timelocked{1}; % JL: ONLY 1 TRIGGER CODE
                save(fullfile(save_path, [params.paradigm '_timelocked']), 'timelocked', '-v7.3');
                clear timelocked

            end
            clear data_ica

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
            data_ica_ds=load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
            params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));

            clear data_ica
            opm_trans = fit_hpi(hpi_path, aux_files{1}, save_path, params);
            % Plot source and head models
            clear headmodels sourcemodel
            headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
            sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;

            timelocked = load(fullfile(save_path, [params.paradigm '_timelocked'])).timelocked;

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

        % Clean up channels

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
                timelocked = load(fullfile(save_path,[params.paradigm '_timelocked.mat'])).timelocked; % JL: LETS TRY THIS
                % mne = fit_mne_opmbci(timelocked,headmodels,sourcemodel,params);
                mne = fit_mne(save_path,timelocked,headmodels,sourcemodel,[],params);    % here i get the filter

                % And now transform activity to source space using the filter
                filt = cellfun(@vecnorm,mne.avg.filter,'UniformOutput',false);
                mne_inv = zeros(length(filt),length(mne.avg.label));
                for i = 1:length(filt)
                    mne_inv(i,:) = filt{i};
                end

                % JL: THERE ARE SEVERAL ISSUES HERE
                % WHY IS mne_inv 1d while the filter is 3D?
                % Am i loading the right data_ica? It has 125 channels when mne_inv seems to have 126. Different subjects?

                % % Plot timelocked
                % h = figure ;
                % plot(timelocked.time,timelocked.avg)
                % title('Motor imagery timelocked')
                % saveas(h, fullfile(save_path, 'figs', 'opm_motorimag_timelocked.jpg'))
                % close all
                % clear mne timelocked
                data_ica = load(fullfile(save_path, [params.paradigm '_data_ica.mat'])).data_ica;
                cfg = [];
                cfg.channel = params.chs;
                data_ica = ft_selectdata(cfg,data_ica);
                data_source = rmfield(data_ica,{'trial','label'});
                for i_trl = 1:length(data_ica.trial)
                    data_source.trial{i_trl} = mapping_matrix*((mne_inv*data_ica.trial{i_trl}));
                end
                data_source.label = sourcemodel.brainstructurelabel;

                 save(fullfile(save_path, [params.paradigm '_data_source']), 'data_source', '-v7.3');
                % clear mne_inv trial time label  data_ica
            end
            clear pos tri roi sourcemodel headmodels mapping_matrix data_ica mne mne_inv timelocked filt
        end 
        close all
    end
end


