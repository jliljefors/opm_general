function fit_mne(save_path, timelocked, headmodel, sourcemodel, sourcemodel_inflated, params)
%fit_mne Performs MNE source reconstruction of data 
%   Expected inputs:
%   save_path: where to save the results
%   timelocked: timelocked/evoked data
%   headmodel: MEG headmodel
%   sourcemodel: sourcemodel
%   sourcemodel_inflated: inflated sourcemodel that can be used for
%                         plotting. If not used set params.plot_inflated 
%                         to false and sourcemodel inflated to []
%   params: parameters struct

if ~isfield(params,'plot_inflated')
    params.plot_inflated = false;
end

%% Prepare leadfields
headmodel = ft_convert_units(headmodel,'cm');
sourcemodel = ft_convert_units(sourcemodel,'cm');
if params.source_fixedori
    sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
end

%% MNE invserse
peak = cell(length(params.trigger_codes),length(params.peaks));

%% Leadfield
cfg = [];
cfg.grad             = timelocked{1}.grad; % sensor positions
cfg.channel          = params.chs;
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
leadfield = ft_prepare_leadfield(cfg,timelocked{1});

for i_trigger = 1:length(params.trigger_codes)
    params.i_trigger = i_trigger;

    %% Set covariance matrix (based on params.noise_cov selection)
    cov = '';
    if isfield(params,'noise_cov')
        if strcmp(params.noise_cov,'all')
            cov = '_covAll';
            timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_all;
        elseif strcmp(params.noise_cov,'resting_state')
            cov = '_covRS';
            if isfield(timelocked{i_trigger},'cov_RS') && size(timelocked{i_trigger}.cov_RS,1) == size(timelocked{i_trigger}.cov,1)
                timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_RS;
            else
                warning('Resting state covariance not existing or incorrect size.');
                break
            end
        elseif strcmp(params.noise_cov,'empty_room')
            cov = '_covER';
            if isfield(timelocked{i_trigger},'cov_ER') && size(timelocked{i_trigger}.cov_ER,1) == size(timelocked{i_trigger}.cov,1)
                timelocked{i_trigger}.cov = timelocked{i_trigger}.cov_ER;
            else
                warning('Empty room covariance not existing or incorrect size.');
                break
            end
        end
    end

    %% Source reconstruction
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = params.chs;
    srcdist = ft_sourceanalysis(cfg, timelocked{i_trigger});
    srcdist.tri = sourcemodel.tri;

    if iscell(srcdist.avg.mom) && size(srcdist.avg.mom{1},1)~=3
        if size(srcdist.avg.mom,1) == 1
            srcdist.avg.mom = cell2mat(srcdist.avg.mom');
        else
            srcdist.avg.mom = cell2mat(srcdist.avg.mom);
        end
    end
    
    h = figure;
    plot(srcdist.time*1e3,srcdist.avg.pow)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.trigger_labels{params.i_trigger}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_trig-' params.trigger_labels{params.i_trigger} '.jpg']))
    close all

    if params.plot_inflated
        srcdist.pos = sourcemodel_inflated.pos;
        srcdist.tri = sourcemodel_inflated.tri;
    end

    for i_peak = 1:length(params.peaks)
        this_peak = FullAreaHalfMax(srcdist,sourcemodel,params.peaks{i_peak}.peak_latency,params);
        this_peak.label = params.peaks{i_peak}.label;

        h = plot_source_distribution(srcdist, this_peak, params);
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' this_peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all

        peak{i_trigger,i_peak} = this_peak;
    end
    clear srcdist

end
save(fullfile(save_path, [params.modality '_mne_peaks' cov]), 'peak'); 

if isfield(params,'save_mne') && params.save_mne
    save(fullfile(save_path, [params.modality '_mne_srcdist' cov]), 'srcdist'); 
end

end