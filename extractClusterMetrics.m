function output_metrics_name = extractClusterMetrics(output_dir,ch_idx_table,file_name,map,Fs)
cd(output_dir)

[~,file_name,~] = fileparts(file_name);
filtered_data_path = fullfile(output_dir,['filt_' file_name '.mda']);
firings_data_path = fullfile(output_dir,['firings_' file_name '.mda']);
% Get cluster metrics from Mountainsort 
MS_metrics_path = fullfile(output_dir,['combine_metrics_new.json']);
fid = fopen(MS_metrics_path); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

m = memmapfile(filtered_data_path,...
    'Offset',20,...
    'Format','int16',...
    'Writable',true);

nCH = size(map,1);
nCOL = numel(m.data)/nCH;


A = readmda(firings_data_path);

primary_channels = A(1,:); 
timestamps = A(2,:)+1; % handle difference between python and matlab
cluster = A(3,:);
cluster_primary_ch = unique([cluster' primary_channels'],'rows');
clear A; 

tBefore=30;
tAfter=45;

valid = timestamps>(tBefore) &  timestamps<(nCOL-tAfter);
timestamps=timestamps(valid);
cluster=cluster(valid);

metrics.date=date;
metrics.nCH = nCH;
metrics.channel_map=map;
metrics.ch_idx_table = ch_idx_table;
clusterID = unique(cluster);
metrics.clusterID=clusterID;
metrics.cluster_info = cell(nCH,1);
metrics.MS_metrics = cell(size(metrics.clusterID));
metrics.waveform=cell(size(metrics.clusterID));
metrics.waveformStd=cell(size(metrics.clusterID));

metrics.P2P = cell(size(metrics.clusterID));
metrics.valley_across_time = cell(size(metrics.clusterID));
metrics.time =cell(size(metrics.clusterID));

metrics.Max_Ins_5s_FR=zeros(size(metrics.clusterID));
metrics.Avg_FR=zeros(size(metrics.clusterID));
metrics.peak_across_time = cell(size(metrics.clusterID));

filt_data = readmda(filtered_data_path);


native_order = ch_idx_table.('Native Intan order');
depth_idx = ch_idx_table.('Depth index');

for clu =1:numel(clusterID)
    
    fprintf('Extracting metrics from cluster #%d...\n',clu)
    metrics.MS_metrics{clu} = val.clusters(clu).metrics;
    primary_ch = cluster_primary_ch(clu,2);
    metrics.cluster_info{clu} = [clu, primary_ch, native_order(primary_ch), depth_idx(primary_ch)];

    p=find(cluster==clusterID(clu));
    unit_firing_times = timestamps(p);
    metrics.time{clu}=timestamps(p);
    
    ISI_in_MS_bins = diff(double(unit_firing_times))/(Fs/1000); % intervals in ms
    [counts] = histcounts(ISI_in_MS_bins,0:1:50);
    [counts_M] = histcounts(ISI_in_MS_bins,0:20:1000);
    [counts_L] = histcounts(ISI_in_MS_bins,0:1000:20000);
    
    second_bins = double(unit_firing_times)/Fs; % what is this
    [counts_5s_1] = histcounts(second_bins,0:5:720);
    [counts_5s_2] = histcounts(second_bins,2.5:5:722.5);
    Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
    
    
    unit=zeros(nCH,tBefore+tAfter+1,length(p),'single');

    for k=1:length(p)
        temp = filt_data(1:nCH,timestamps(p(k))-tBefore:timestamps(p(k))+tAfter);
        unit(:,:,k)=1/10*single(temp);
    end
    % valley_across_time = min(reshape(unit,[size(unit,1)*size(unit,2),size(unit,3)]));
    
    mean_plot=mean(unit,3);
    if size(unit,3)>3.6E6
        STD_plot=std(unit(:,:,1:1E5),0,3);
    else
        STD_plot=std(unit,0,3);
    end
    
    P2P = max(mean_plot')-min(mean_plot');
    
    ch_with_max_P2P=find(P2P==max(P2P));
    ch_with_max_P2P=ch_with_max_P2P(1);
    
    [sortP2P,sP2P]=sort(P2P,'descend');
    peak_across_time = unit(primary_ch,tBefore+1,:);
    valley_across_time = unit(ch_with_max_P2P,tBefore+1,:);
    valley_across_time  =valley_across_time(:);
    peak_across_time =peak_across_time(:);
    
    metrics.waveform{clu}=mean_plot;
    metrics.P2P{clu}=P2P;
    metrics.waveformStd{clu}=STD_plot;
    maxV = max(abs(peak_across_time),abs(peak_across_time));
    if maxV<3276
        metrics.valley_across_time{clu}=int16(valley_across_time);
        metrics.peak_across_time{clu}=int16(peak_across_time);
    else
        metrics.valley_across_time{clu}=valley_across_time;
        metrics.peak_across_time{clu}=peak_across_time;
    end
    
    metrics.Max_Ins_5s_FR(clu)= Max_Ins_5s_FR;
    metrics.Avg_FR(clu)= sum(cluster==clusterID(clu))/(nCOL/20E3);

    mean_plot=mean(unit,3);
    if size(unit,3)>3.6E6
        STD_plot=std(unit(:,:,1:1E5),0,3);
    else
        STD_plot=std(unit,0,3);
    end
    P2P = max(mean_plot')-min(mean_plot');
    ch_with_max_P2P=find(P2P==max(P2P));
    ch_with_max_P2P=ch_with_max_P2P(1);
    [sortP2P,sP2P]=sort(P2P,'descend');
    peak_across_time = unit(primary_ch,tBefore+1,:);
    valley_across_time = unit(ch_with_max_P2P,tBefore+1,:);
    valley_across_time  =valley_across_time(:);
    peak_across_time    =peak_across_time(:);
    

end
metrics.list=clusterID;
metrics.cluster_primary_ch=cluster_primary_ch;
fprintf('Cluster metric extraction complete!\n')

output_metrics_name = fullfile(output_dir,[metrics.date '_cluster_metrics']);

save(output_metrics_name,'metrics')
end