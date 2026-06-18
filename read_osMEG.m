function [data_cleaned] = read_osMEG(opm_file, aux_file, save_path, params, trialinfo)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% opm_file: path to the file containing the OPM-MEG data
% aux_file: path to synchronously recorded file containing ECG, EOG and/or
%           EEG. If not appliccable set to []
% save_path: path where to save results
% params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

if ~exist(opm_file,'file')
    error(['Did not find OPM file: ' opm_file])
end

if isempty(aux_file)
    opm_only = true;
else
    opm_only = false;
    if ~exist(aux_file,'file')
        error(['Did not find AUX file: ' aux_file])
    end
end

% get opm files
[folder, base, ext] = fileparts(opm_file);
files = dir(fullfile(folder, [base '*' ext]));
names = {files.name};
splitIdx = zeros(numel(names),1);  % original file = 0
for k = 1:numel(names)
    tok = regexp(names{k}, [regexptranslate('escape', base) '-(\d+)' regexptranslate('escape', ext) '$'], ...
        'tokens', 'once');
    if ~isempty(tok)
        splitIdx(k) = str2double(tok{1});
    end
end
[~, order] = sort(splitIdx);
files = files(order);
opm_files = fullfile({files.folder}, {files.name})';

%% --- Read triggers ---
% OPM
trl_opm=[];
cfg = [];
% cfg.datafile        = opm_file;
cfg.dataset = opm_files;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);
i_trig_opm = find(contains(opm_raw.label,'di'));
trig = opm_raw.trial{1}(i_trig_opm,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_opm(:,1) = find(trig)-(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,2) = find(trig)+(params.post+params.pad)*opm_raw.fsample;
trl_opm(:,3) = -(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,4) = opm_raw.trial{1}(i_trig_opm,trig);
trl_opm(:,1:2) = trl_opm(:,1:2) + floor(params.delay*opm_raw.fsample); % adjust for stim delay
trl_opm = round(trl_opm);

trl_opm = trl_opm(ismember(trl_opm(:,4), unique(cell2mat(params.trigger_codes))),:);

if ~opm_only
    % AUX
    trl_aux=[];
    cfg = [];
    cfg.datafile        = aux_file;
    aux_raw = ft_preprocessing(cfg);
    i_trig_aux = find(contains(aux_raw.label,'STI101'));
    trig = aux_raw.trial{1}(i_trig_aux,:)>0.5;
    trig = [false trig(2:end)&~trig(1:end-1)];
    trl_aux(:,1) = find(trig)-(params.pre+params.pad)*aux_raw.fsample;
    trl_aux(:,2) = find(trig)+(params.post+params.pad)*aux_raw.fsample;
    trl_aux(:,3) = -(params.pre+params.pad)*aux_raw.fsample;
    trl_aux(:,4) = aux_raw.trial{1}(i_trig_aux,trig);
    trl_aux(:,1:2) = trl_aux(:,1:2) + floor(params.delay*aux_raw.fsample); % adjust for stim delay
    trl_aux = round(trl_aux);
    
    % Check if uneven amount of trial. If so assume error in beginning.
    if size(trl_aux,1) > size(trl_opm,1)
        trl_aux = trl_aux((end-size(trl_opm,1)+1):end,:);
    elseif size(trl_aux,1) < size(trl_opm,1)
        trl_opm = trl_opm((end-size(trl_aux,1)+1):end,:);
    end
    if trl_aux(:,4) ~= trl_opm(:,4) % Throw error if trials don't match.
        error('events do not match')
    end
    
    %% AUX data filter & epoch
    cfg = [];
    if isfield(params.filter,'hp_freq') && ~isempty(params.filter.hp_freq)
        cfg.hpfilter        = 'yes';         
        cfg.hpfreq          = params.filter.hp_freq;
        cfg.hpinstabilityfix  = 'reduce';
    end
    if isfield(params.filter,'lp_freq') && ~isempty(params.filter.lp_freq)
        cfg.lpfilter        = 'yes';         
        cfg.lpfreq          = params.filter.lp_freq;
    end
    if isfield(params.filter,'bp_freq') && ~isempty(params.filter.bp_freq)
        cfg.bpfilter        = 'yes';         
        cfg.bpfreq          = params.filter.bp_freq;
        cfg.bpinstabilityfix  = 'reduce';
    end
    aux_epo = ft_preprocessing(cfg,aux_raw);
    
    cfg = [];
    cfg.trl             = trl_aux;
    aux_epo = ft_redefinetrial(cfg,aux_epo);
    
    cfg = [];
    cfg.dftfilter       = 'yes';        
    cfg.dftfreq         = params.filter.notch;
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-params.pre 0];
    aux_epo = ft_preprocessing(cfg,aux_epo);
elseif any(contains(opm_raw.label,'ai173'))
    ecg_ch = find(contains(opm_raw.label,'ai173'));
    eog1_ch = find(contains(opm_raw.label,'ai174'));
    eog2_ch = find(contains(opm_raw.label,'ai175'));
    disp('Changing channel types for ECG/EOG')
    if ~isempty(ecg_ch)
        opm_raw.label{ecg_ch} = 'ECG';
        opm_raw.hdr.label{ecg_ch} = 'ECG';
        opm_raw.hdr.chantype{ecg_ch} = 'ecg';
    end
    if ~isempty(eog1_ch)
        opm_raw.label{eog1_ch} = 'EOG001';
        opm_raw.hdr.label{eog1_ch} = 'EOG1';
        opm_raw.hdr.chantype{eog1_ch} = 'eog';
    end
    if ~isempty(eog2_ch)
        opm_raw.label{eog2_ch} = 'EOG002';
        opm_raw.hdr.label{eog2_ch} = 'EOG2';
        opm_raw.hdr.chantype{eog2_ch} = 'eog';
    end
end

%% OPM data filter & epoch
cfg = [];
if isfield(params.filter,'hp_freq') && ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes';         
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
end
if isfield(params.filter,'lp_freq') && ~isempty(params.filter.lp_freq)
    cfg.lpfilter        = 'yes';         
    cfg.lpfreq          = params.filter.lp_freq;
end
if isfield(params.filter,'bp_freq') && ~isempty(params.filter.bp_freq)
    cfg.bpfilter        = 'yes';         
    cfg.bpfreq          = params.filter.bp_freq;
    cfg.bpinstabilityfix  = 'reduce';
end
opm_epo = ft_preprocessing(cfg, opm_raw);

cfg = [];
cfg.trl             = trl_opm;
opm_epo = ft_redefinetrial(cfg,opm_epo);

cfg = [];
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_epo = ft_preprocessing(cfg, opm_epo);

% Resample
if ~opm_only
    cfg            = [];
    cfg.time = aux_epo.time;
    cfg.detrend    = 'no';
    opm_epo_ds = ft_resampledata(cfg, opm_epo);
else
    cfg            = [];
    cfg.resamplefs = 1000;
    cfg.detrend    = 'no';
    opm_epo_ds = ft_resampledata(cfg, opm_epo);
end

%% Combine data
if ~opm_only
    EOG_channels = find(contains(aux_epo.label,'EOG'));
    ECG_channels = find(contains(aux_epo.label,'ECG'));
    EEG_channels = find(contains(aux_epo.label,'EEG'));
    MISC_channels = find(contains(aux_epo.label,'MISC'));
    TRIG_channels = find(contains(aux_epo.label,'STI101'));
    include_channels = [EOG_channels; ECG_channels; EEG_channels; MISC_channels; TRIG_channels];
    
    data = opm_epo_ds; 
    data.elec = aux_epo.elec;
    data.time = aux_epo.time;
    data.label = [data.label; aux_epo.label(include_channels)];
    data.hdr.label = data.label;
    data.hdr.nChans = data.hdr.nChans + length(include_channels);
    data.hdr.chantype = [data.hdr.chantype; aux_epo.hdr.chantype(include_channels)];
    data.hdr.chanunit = [data.hdr.chanunit; aux_epo.hdr.chanunit(include_channels)];
    data.sampleinfo = aux_epo.sampleinfo;
    data.trialinfo = aux_epo.trialinfo;
    n_smpl = size(aux_epo.trial{1},2);
    for i = 1:length(data.trial)
        data.trial{i} = [data.trial{i}(:,1:n_smpl); aux_epo.trial{i}(include_channels,:)]; 
    end
else
    data = opm_epo_ds;
end

data.trialinfo = trialinfo;
data.trialinfo.orig_trial_number = (1:height(data.trialinfo))';
n_trials_initial = height(data.trialinfo);

%% Find & remove bad opm channels
[badchs, badchs_flat, badchs_std, badchs_neighbors, badchs_outlier] = opm_badchannels(opm_raw, trl_opm, params);
save(fullfile(save_path, [params.paradigm '_badchs']), ...
    'badchs_flat', ...
    'badchs_std', ...
    'badchs_neighbors', ...
    'badchs_outlier',"-v7.3"); 

cfg = [];
cfg.channel = setdiff(data.label,badchs);
data = ft_selectdata(cfg, data);
data.coordsys = 'scanras';
data.grad.coordsys = 'scanras';

%% Spatiotemporal filtering
cfg = []; % separate ExG channels
cfg.channel = {'EOG*', 'ECG*'};
ExG = ft_selectdata(cfg,data);

if params.do_hfc
    cfg = [];
    cfg.channel = '*bz';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    data_cleaned = ft_denoise_hfc(cfg, data);
elseif params.do_amm
    cfg = [];
    cfg.channel = '*bz';
    cfg.updatesens = 'yes';
    cfg.residualcheck = 'no';
    cfg.amm = [];
    cfg.amm.order_in = params.amm_in;
    cfg.amm.order_out = params.amm_out;
    cfg.amm.thr = params.amm_thr;
    cfg.trials = true(1,numel(data.trial)); cfg.trials(end) = false;
    data_cleaned = ft_denoise_amm(cfg, data);
else
    data_cleaned = data;
end

% Recombine with ExG channels
data_cleaned.label = vertcat(data_cleaned.label,ExG.label);
data_cleaned.hdr = data.hdr;
incl = ismember(data.hdr.label,data_cleaned.label);
data_cleaned.hdr.label = data.hdr.label(incl);
data_cleaned.hdr.chantype = data.hdr.chantype(incl);
data_cleaned.hdr.chanunit = data.hdr.chanunit (incl);
for i = 1:length(data_cleaned.trial)
    data_cleaned.trial{i} = vertcat(data_cleaned.trial{i}, ExG.trial{i}); 
end

%% Visual trial rejection
cfg         = [];
cfg.method  = 'summary';
cfg.channel = {'*bz'};
cfg.metric  = 'std';
cfg.layout  = params.layout;
data_cleaned = ft_rejectvisual(cfg, data_cleaned);

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) 
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    data_cleaned = ft_resampledata(cfg, data_cleaned);
end

%% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
data_cleaned = ft_selectdata(cfg, data_cleaned); 

%% Convert to sensor definitions to cm
data_cleaned.grad = ft_convert_units(data_cleaned.grad,'cm');

%% Preprocessing info (carried through pipeline; consolidated in subject_info.mat)
data_cleaned.preproc_info = struct();
data_cleaned.preproc_info.badchs           = badchs;
data_cleaned.preproc_info.badchs_flat      = badchs_flat;
data_cleaned.preproc_info.badchs_std       = badchs_std;
data_cleaned.preproc_info.badchs_neighbors = badchs_neighbors;
data_cleaned.preproc_info.badchs_outlier   = badchs_outlier;
data_cleaned.preproc_info.n_trials_initial = n_trials_initial;
if istable(data_cleaned.trialinfo) && any(strcmp(data_cleaned.trialinfo.Properties.VariableNames, 'orig_trial_number'))
    data_cleaned.preproc_info.bad_trials = setdiff((1:n_trials_initial)', data_cleaned.trialinfo.orig_trial_number);
else
    data_cleaned.preproc_info.bad_trials = [];
    warning('read_osMEG: orig_trial_number column not found on data_cleaned.trialinfo; bad_trials not tracked.');
end

% %% Save bad trials JL DOESNT WORK
% [~,idx]=ismember(data_cleaned.sampleinfo,badtrl_jump,'rows');
% badtrl_opm_jump = find(idx);
% [~,idx]=ismember(data_cleaned.sampleinfo,badtrl_std,'rows');
% badtrl_opm_std = find(idx);
% save(fullfile(save_path, [params.sub '_opm_badtrls']), ...
%     'badtrl_opm_jump', ...
%     'badtrl_opm_std',"-v7.3"); 

end