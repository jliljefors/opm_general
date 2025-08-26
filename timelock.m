function [timelocked] = timelock(data, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
amp_label     = ft_getopt(params, 'amp_label', 'na');
amp_scaler     = ft_getopt(params, 'amp_scaler', 1);

timelocked = cell(length(params.trigger_codes),1);

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
    % Average trials
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = [-params.pre 0];
    % cfg.trials = find(data.trialinfo==params.trigger_codes(i_trigger));      %  JL: This clashes with my trialinfo, and is irrelevant since i only have one paradigm
    cfg.channel = params.include_chs;                                        %  JL: Exclude bad channels!
    timelocked{i_trigger} = ft_timelockanalysis(cfg, data);
    timelocked{i_trigger}.trigger_code = params.trigger_codes(i_trigger);
    timelocked{i_trigger}.trigger_label = params.trigger_labels(i_trigger);

    % Butterfly plot
    chs = find(contains(timelocked{i_trigger}.label,ft_channelselection(params.chs,timelocked{i_trigger}.label)));
    h = figure;
    plot(timelocked{i_trigger}.time*1e3,timelocked{i_trigger}.avg(chs,:)*amp_scaler)
    xlabel('t [msec]')
    ylabel(amp_label)
    %xlim([-params.pre params.post]*1e3);
    title(['Evoked - ' params.trigger_labels{i_trigger} ' (n_{trls}=' num2str(length(timelocked{i_trigger}.cfg.trials)) ')'])
    saveas(h, fullfile(save_path, 'figs', [params.paradigm '_butterflyPlot_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
    close all

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
        %xlim([-params.pre params.post]*1e3);
        saveas(h, fullfile(save_path, 'figs', [params.paradigm '_evoked-' params.plot_channel '_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
        close all
    end

    if isfield(params,'plot_latency') % latency to plot topogrpaphy at in seconds
        cfg = [];
        cfg.xlim = [params.plot_latency-0.010 params.plot_latency+0.010];
        cfg.layout = params.layout; 
        cfg.parameter = 'avg';
        h = figure;
        ft_topoplotER(cfg, timelocked{i_trigger});
        axis on
        colorbar
        saveas(h, fullfile(save_path, 'figs', [params.paradim '_topography_trig-' params.trigger_labels{i_trigger} '_' params.modality '.jpg']))
        close all
    end
end

end