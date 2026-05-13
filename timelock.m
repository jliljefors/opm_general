function [timelocked] = timelock(data, save_path, params)
%Timelock Average trial/epoch data and plot
%   Averages trials/epochs for each trigger defined in params.trigger_codes
%   and plots a butterfly (with global field power). If a plot_channel or
%   plot_latency are selected in params the function will generate a plot
%   of the channel with variation over trials and topography of the
%   latency, respectively. For the topogrpahy a layout should also be
%   defined in params. 
%   Required inputs:
%   data: data containg individual trials
%   save_path: path where to save the generated figures
%   params: parameter struct

amp_label     = ft_getopt(params, 'amp_label', 'na'); % used for yaxis of the figures generated (ex: 'B [fT]')
amp_scaler     = ft_getopt(params, 'amp_scaler', 1); % scaling factor used for figures (ex: 1e15 to display signals in fT)

timelocked = cell(length(params.trigger_codes),1);

if isfield(params,'trigger_codes') && isempty(params.trigger_codes)
    params.trigger_codes = num2cell(unique(data.trialinfo));
    for i = 1:length(params.trigger_codes)
        params.trigger_labels{i} = num2str(params.trigger_codes{i});
    end
end

% Cut off padding
cfg = [];
cfg.channel = params.chs;
cfg.latency = [-params.pre params.post];
data = ft_selectdata(cfg, data);

% Demean
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-params.pre 0];
data = ft_preprocessing(cfg,data);

for i_trigger = 1:length(params.trigger_codes)
    % Select trials
    if numel(params.trigger_codes) == 1
        trls = true(1,height(data.trialinfo));
    elseif isnumeric(params.trigger_codes{i_trigger}) && length(params.trigger_codes{i_trigger})==1 % trigger code
        trls = find(data.trialinfo==params.trigger_codes{i_trigger});
    elseif isnumeric(params.trigger_codes{i_trigger}) && length(params.trigger_codes{i_trigger})>1 % list of trials
        trls = params.trigger_codes{i_trigger};
    elseif ischar(params.trigger_codes{i_trigger}) && strcmp(params.trigger_codes{i_trigger},'all') %
        trls = 1:length(data.trial);
    end

    % Average trials
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = [-params.pre 0];
    cfg.trials = trls;
    timelocked{i_trigger} = ft_timelockanalysis(cfg, data);
    timelocked{i_trigger}.trigger_code = params.trigger_codes(i_trigger);
    timelocked{i_trigger}.trigger_label = params.trigger_labels(i_trigger);

    %%
    gfp = std(timelocked{i_trigger}.avg,0,1)*params.amp_scaler;
    
    %% Plot butterfly
    h = figure;
    left = 0.1;
    bottom = 0.1;
    width = 0.8;
    gap = 0.10;  % gap between plots
    height_total = 0.8 - gap;  
    height_top = 3/4 * height_total;
    height_bottom = 1/4 * height_total;
    
    % Butterfly
    ax1 = axes('Position', [left, bottom + height_bottom + gap, width, height_top]);
    plot(timelocked{i_trigger}.time*1e3,timelocked{i_trigger}.avg*params.amp_scaler)
    xlabel('t [msec]')
    ylabel(params.amp_label)
    xlim([-params.pre params.post]*1e3);
    title(['Evoked ' params.modality ' - ' params.trigger_labels{i_trigger} ' (n_{trls}=' num2str(length(timelocked{i_trigger}.cfg.trials)) ')'])

    % GFP
    ax2 = axes('Position', [left, bottom, width, height_bottom]);
    plot(timelocked{i_trigger}.time*1e3,gfp,'k')
    ax2.XTickLabel = [];
    ylabel('GFP')
    xlim([-params.pre params.post]*1e3);
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all

    %% Plot selected channel
    if isfield(params,'plot_channel') && sum(contains(timelocked{i_trigger}.label,params.plot_channel)) == 1 % only if a single channel is selected
        i_plot_ch = find(contains(timelocked{i_trigger}.label,params.plot_channel)); % pick like 'L204' or 'L204_bz' 
        h = figure;
        hold on
        for i_trl = find(data.trialinfo==params.trigger_codes(i_trigger))'
            plot(data.time{i_trl}*1e3, data.trial{i_trl}(i_plot_ch,:)*params.amp_scaler,'Color',[211 211 211]/255)
        end
        plot(timelocked{i_trigger}.time*1e3, timelocked{i_trigger}.avg(i_plot_ch,:)*amp_scaler,'Color',[0 0 0]/255)
        hold off
        title(['Channel: ' params.plot_channel])
        ylabel(amp_label)
        xlabel('time [ms]')
        xlim([-params.pre params.post]*1e3);
        saveas(h, fullfile(save_path, 'figs', [params.paradigm '_evoked-' params.plot_channel '_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
        close all
    end

    %% Plot topogrpaphy at selected time
    if isfield(params,'plot_latency') % latency to plot topogrpaphy at in seconds
        cfg = [];
        cfg.xlim = [params.plot_latency-0.010 params.plot_latency+0.010];
        cfg.layout = params.layout; 
        cfg.parameter = 'avg';
        cfg.comment = 'no';
        h = figure;
        ft_topoplotER(cfg, timelocked{i_trigger});
        axis on
        colorbar
        saveas(h, fullfile(save_path, 'figs', [params.paradim '_topography_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
        close all
    end
end

end