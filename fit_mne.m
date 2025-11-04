function fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelocked, headmodels, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(params,'plot_inflated') && params.plot_inflated && ~isempty(sourcemodel_inflated)
    params.plot_inflated = false;
end

%% Prepare leadfields
headmodel = headmodels.headmodel_meg;
headmodel = ft_convert_units(headmodel,'cm');
sourcemodel = ft_convert_units(sourcemodel,'cm');
sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';

%% MNE invserse
squidmag_peak = cell(length(params.trigger_code),length(params.peaks));
squidgrad_peak = cell(length(params.trigger_code),length(params.peaks));
opm_peak = cell(length(params.trigger_code),length(params.peaks));

%% Leadfields
cfg = [];
cfg.grad             = squidmag_timelocked{1}.grad; % sensor positions
cfg.channel          = 'megmag';
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
cfg.normalize        = 'yes';
leadfield_squidmag = ft_prepare_leadfield(cfg,squidmag_timelocked{1});

cfg = [];
cfg.grad             = squidgrad_timelocked{1}.grad; % sensor positions
cfg.channel          = 'megplanar';
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
cfg.normalize        = 'yes';
leadfield_squidgrad = ft_prepare_leadfield(cfg,squidgrad_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad; % sensor positions
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
cfg.normalize        = 'yes';
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

for i_trigger = 1:length(params.trigger_code)
    
    %% MEG-MAG
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_squidmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'megmag';
    mne_squidmag = ft_sourceanalysis(cfg, squidmag_timelocked{i_trigger});
    mne_squidmag.tri = sourcemodel.tri;

    params.modality = 'squidmag';

    if iscell(mne_squidmag.avg.mom)
        if size(mne_squidmag.avg.mom{1},1)>1
            for i = 1:length(mne_squidmag.avg.mom)
                mne_squidmag.avg.mom{i} = vecnorm(mne_squidmag.avg.mom{i},1);
            end
        elseif size(mne_squidmag.avg.mom,1) == 1
            mne_squidmag.avg.mom = cell2mat(mne_squidmag.avg.mom');
        else
            mne_squidmag.avg.mom = cell2mat(mne_squidmag.avg.mom);
        end
    end
    h = figure;
    plot(mne_squidmag.time*1e3,mne_squidmag.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_mne_sourcepow_trig-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        squidmag_peak{i_trigger,i_peak} = FullAreaHalfMax(mne_squidmag,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        mne_squidmag.pos = sourcemodel_inflated.pos;
        mne_squidmag.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidmag_peak{i_trigger,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, mne_squidmag)
        lighting gouraud
        material dull
        title(['SQUID-MAG (FAHM=' num2str(squidmag_peak{i_trigger,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidmag_peak{i_trigger,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_' squidmag_peak{i_trigger,i_peak}.label '_mne_trig' params.phalange_labels{i_trigger} '.jpg']))
        close all
    end

    %% MEG-GRAD
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield_squidgrad;
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'megplanar';
    mne_squidgrad = ft_sourceanalysis(cfg, squidgrad_timelocked{i_trigger});
    mne_squidgrad.tri = sourcemodel.tri;

    params.modality = 'squidgrad';

    if iscell(mne_squidgrad.avg.mom)
        if size(mne_squidgrad.avg.mom{1},1)>1
            for i = 1:length(mne_squidgrad.avg.mom)
                mne_squidgrad.avg.mom{i} = vecnorm(mne_squidgrad.avg.mom{i},1);
            end
        elseif size(mne_squidgrad.avg.mom,1) == 1
            mne_squidgrad.avg.mom = cell2mat(mne_squidgrad.avg.mom');
        else
            mne_squidgrad.avg.mom = cell2mat(mne_squidgrad.avg.mom);
        end
    end
    h = figure;
    plot(mne_squidgrad.time*1e3,mne_squidgrad.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_mne_sourcepow_trig-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        squidgrad_peak{i_trigger,i_peak} = FullAreaHalfMax(mne_squidgrad,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        mne_squidgrad.pos = sourcemodel_inflated.pos;
        mne_squidgrad.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidgrad_peak{i_trigger,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, mne_squidgrad)
        lighting gouraud
        material dull
        title(['SQUID-GRAD (FAHM=' num2str(squidgrad_peak{i_trigger,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidgrad_peak{i_trigger,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_' squidgrad_peak{i_trigger,i_peak}.label '_mne_trig' params.phalange_labels{i_trigger} '.jpg']))
        close all
    end

    %% OPM
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = '*bz';
    mne_opm = ft_sourceanalysis(cfg, opm_timelocked{i_trigger});
    mne_opm.tri = sourcemodel.tri;

    params.modality = 'opm';

    if iscell(mne_opm.avg.mom)
        if size(mne_opm.avg.mom{1},1)>1
            for i = 1:length(mne_opm.avg.mom)
                mne_opm.avg.mom{i} = vecnorm(mne_opm.avg.mom{i},1);
            end
        elseif size(mne_opm.avg.mom,1) == 1
            mne_opm.avg.mom = cell2mat(mne_opm.avg.mom');
        else
            mne_opm.avg.mom = cell2mat(mne_opm.avg.mom);
        end
    end
    h = figure;
    plot(mne_opm.time*1e3,mne_opm.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_mne_sourcepow_trig-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all
    
    for i_peak = 1:length(params.peaks)
        opm_peak{i_trigger,i_peak} = FullAreaHalfMax(mne_opm,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        mne_opm.pos = sourcemodel_inflated.pos;
        mne_opm.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)  
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = opm_peak{i_trigger,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, mne_opm)
        lighting gouraud
        material dull
        title(['OPM (FAHM=' num2str(opm_peak{i_trigger,i_peak}.fahm,3) 'cm^2; t=' num2str(round(opm_peak{i_trigger,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_opm_' opm_peak{i_trigger,i_peak}.label '_mne_trig' params.phalange_labels{i_trigger} '.jpg']))
        close all
    end

end

end