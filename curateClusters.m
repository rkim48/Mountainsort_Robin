% Curation of clusters with user defined parameters
metrics = load([output_metrics_name '.mat']);
metrics =metrics.metrics;

for unit=1:numel(metrics.list)
    
    wav = metrics.waveform{unit};
    wav_std = metrics.waveformStd{unit};

    P2P = min(wav');
    ch=find(P2P==min(P2P));
    ch=ch(1);
    tp=find(wav(ch,:)==min(P2P));
    tp=tp(1);
    metrics.peakCV(unit)=abs(wav_std(ch,tp)/min(P2P));

    unit_firing_times = metrics.time{unit};
    ISI_in_MS_bins = diff(double(unit_firing_times))/(Fs/1E3); % Intervals in miliseconds. 
    [counts] = histcounts(ISI_in_MS_bins,0:2.5:200); % 0 to 200 ms with 2.5 ms bins (paper below)
    % https://www.nature.com/articles/s41598-019-38924-w#Fig2

    metrics.firings_first_bin = counts(1);
    metrics.single(unit) = counts(1)/numel(unit_firing_times)<=0.01;
    metrics.violation(unit) = counts(1)/numel(unit_firing_times);
    
    max_time = max(cell2mat(metrics.time))/Fs;
    metrics.avg_FR(unit)=numel(unit_firing_times)/max_time;

% we allow one spike per min.
end

output_curated_name = fullfile(output_dir,[metrics.date '_curated_cluster_metrics']);
save(output_curated_name,'metrics')