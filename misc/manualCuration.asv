% Visualizing curated clusters
function visualizeClusters(output_curated_name,output_dir,Fs,animalID)

metrics = load([output_curated_name '.mat']);
metrics =metrics.metrics;
date_str = metrics.date;

nCH = metrics.nCH;
mkdir(fullfile(output_dir,'figures'))
accept_dir = fullfile(output_dir,'accepted');
reject_dir = fullfile(output_dir,'rejected');
idk_dir = fullfile(output_dir,'idk');
if ~isfolder(accept_dir)
    mkdir(accept_dir);
end
if ~isfolder(reject_dir)
    mkdir(reject_dir);
end
if ~isfolder(idk_dir)
    mkdir(idk_dir);
end

for unit = 1:nCH
    figure('WindowState','maximized');
    pause(1);
    fig_rows = 3;
    fig_cols = 4;

    
    cluster_info = metrics.cluster_info{unit};
   
    % ISI histogram 
    unit_firing_times = metrics.time{unit};
    amplitudes_at_firing_times = metrics.peak_across_time{unit};
    ISI_in_MS_bins = diff(double(unit_firing_times))/(Fs/1E3); % Intervals in miliseconds. 
    [counts] = histcounts(ISI_in_MS_bins,0:2.5:200); 
    subplot(fig_rows,fig_cols,2:fig_cols);
    histogram('BinEdges',0:2.5:200,'BinCounts',counts);
    xlabel('Time (ms)')
    ylabel('Frequency');
    title('ISI distribution');
    
    % Templates along shank  
    plot_color = [0 0.4470 0.7410];
    subplot(fig_rows,fig_cols,1:fig_cols:fig_cols*fig_rows)
    min_val = 0; max_val = 0;
    for ch = 1:nCH
        offset = metrics.ch_idx_table{ch,2};
        template = metrics.waveform{unit}(ch,:);
        X = 1:size(template,2);
        template_w_offset = template + (32 - offset)*50; % hard coded
        if min(template_w_offset) < min_val
            min_val = min(template_w_offset);
        end
        if max(template_w_offset) > max_val
            max_val = max(template_w_offset);
        end
        plot(X,template_w_offset,'color',plot_color)
        hold on;
    end
    ylim([min_val-70, max_val+20])
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    obj = scalebar('XLen',20,'XUnit','ms','YUnit','uV'); 
    obj.Position = [5,min_val-30]; 
    obj.hTextX_Pos = [2,-15];
    obj.hTextY_Pos = [-2,0];
    title('Templates depth-ordered');
    
    % Amplitude vs. time
    subplot(fig_rows,fig_cols,fig_cols+(2:fig_cols));
    scatter(unit_firing_times/Fs, amplitudes_at_firing_times,10,'filled');
    xlabel('Time (s)')
    ylabel('Peak amplitude')
    title('Event distribution')
    
    % Other notable metrics
    avg_FR = metrics.avg_FR(unit);
    violations = metrics.violation(unit);
    unit_type = metrics.unit_type{unit};
    title_str = 'Cluster metrics';
%     metric_keys = {'Firing rate','Refractory violations 2.5ms','Unit type'};
%     metric_values = {avg_FR, violations, unit_type};
    subplot(fig_rows,fig_cols,(fig_rows-1)*fig_cols+2);
    y_start = 0.9; x_start = 0.1; line_spacing = 0.2; font_size = 12;
    
    text(x_start,y_start, title_str,'FontWeight','Bold','FontSize',font_size); 
    text(x_start,y_start - line_spacing, sprintf('Firing rate: %.2f Hz',avg_FR),'FontSize',font_size); 
    text(x_start,y_start - 2 * line_spacing, sprintf('Refractory violations 2.5ms: %.2f%%',violations*100),'FontSize',font_size);
    text(x_start,y_start - 3 * line_spacing, sprintf('Unit type: %s',unit_type),'FontSize',font_size); 

    axis off
    
    % Overlapping templates 
    subplot(fig_rows,fig_cols,(fig_rows-1)*fig_cols+(fig_cols-1:fig_cols));
    plot(X,metrics.waveform{unit}')
    ylabel('Amplitude (uV)')
    xlabel('Time (ms)')
    title('Overlapping templates')
    
   if strcmp(unit_type,'Single unit')
        single_or_multi= 'Single ';
    else
        single_or_multi = 'Multi-';
   end
    
    title_str = sprintf('%s %s: Unit %d',date_str,animalID,unit);
    sgtitle(title_str,'fontweight','bold');
    image = fullfile(output_dir,'figures',sprintf('%sunit_%d.png',single_or_multi,unit));
    saveas(gcf,image)
    
    while 1
        fprintf('Press y/n/d for accept/reject/idk\n');       
        w=waitforbuttonpress;
        if w
            value = get(gcf, 'CurrentCharacter');
        end
        switch value
            case 121         
                fprintf('Yes\n');
                text(0.5,0.5,'yes','Color','g','Fontsize',20);
                movefile(image,accept_dir)
                break
            case 110
                fprintf('No\n');
                text(0.5,0.5,'No','Color','r','Fontsize',20);
                movefile(image,reject_dir)
                break
            case 100
                fprintf("Don't know\n");
                text(0.5,0.5,"Don't know",'Color','blue','Fontsize',20);
                movefile(image,idk_dir)
                break
            otherwise
                fprintf('Error: Press y/n\n');
        end
        pause(0.5);
        close;
    end
    
end
fprintf('Cluster visualization complete!\n')
end