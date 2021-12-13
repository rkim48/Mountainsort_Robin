% Visualizing curated clusters

metrics = load([output_curated_name '.mat']);
metrics =metrics.metrics;
mkdir(fullfile(output_dir,'figures'))
for unit = 1:numel(metrics.list)
    figure;
    subplot(2,1,1);
    if metrics.single(unit)
        single_or_multi = 'Single ';
    else
        single_or_multi = 'Multi-';
    end
    plot(metrics.waveform{unit}')
    
    unit_firing_times = metrics.time{unit};
    ISI_in_MS_bins = diff(double(unit_firing_times))/(Fs/1E3); % Intervals in miliseconds. 
    [counts] = histcounts(ISI_in_MS_bins,0:2.5:200); 
    subplot(2,1,2);
    histogram('BinEdges',0:2.5:200,'BinCounts',counts);
    avg_FR = metrics.avg_FR(unit);
    
    title_str = sprintf('%sunit %d with average firing rate: %.1f spikes/s',single_or_multi,unit,avg_FR);
    sgtitle(title_str);
    
    saveas(gcf,fullfile(output_dir,'figures',sprintf('%sunit_%d.png',single_or_multi,unit)))
    close
end