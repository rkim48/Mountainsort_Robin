% Curation of clusters with user defined parameters


cluster_metrics = load([output_metrics_name '.mat']);

cluster_metrics =cluster_metrics.metrics;

for unit=1:numel(cluster_metrics.list)
    
    wav = cluster_metrics.waveform{unit};
    wav_std = cluster_metrics.waveformStd{unit};

    P2P = min(wav');
    ch=find(P2P==min(P2P));
    ch=ch(1);
    tp=find(wav(ch,:)==min(P2P));
    tp=tp(1);
    cluster_metrics.peakCV(unit)=abs(wav_std(ch,tp)/min(P2P));

    FiringTimeForThisUnit = cluster_metrics.time{unit};
    ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/(Fs/1E3); % Intervals in miliseconds. 
    [counts] = histcounts(ISI_in_MS_bins,0:2.5:200); % 0 to 200 ms with 2.5 ms bins (paper below)
    % https://www.nature.com/articles/s41598-019-38924-w#Fig2
    maxTime = max(cell2mat(cluster_metrics.time))/Fs;
   
    Sec_bins = double(FiringTimeForThisUnit)/Fs;
    [counts_5s_1] = histcounts(Sec_bins,0:5:720);
    [counts_5s_2] = histcounts(Sec_bins,2.5:5:722.5);
    Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;

    cluster_metrics.Max_Ins_5s_FR(unit)=Max_Ins_5s_FR;
    cluster_metrics.Avg_FR(unit)=numel(FiringTimeForThisUnit)/maxTime;

    try
        cluster_metrics.multi1(unit)= counts(1);
        cluster_metrics.multi2(unit)= counts(2);
        cluster_metrics.single(unit)= ((counts(2)+counts(1))/numel(FiringTimeForThisUnit))<=0.01;
        cluster_metrics.violation(unit)= ((counts(2)+counts(1))/numel(FiringTimeForThisUnit));
    catch
        try
            cluster_metrics.multi2(unit)= counts(2);
        catch
        end
    end
% we allow one spike per min.
end
result = cluster_metrics;
result.param=param;
result.white = white;
result.msp=msp;
result.sftest=sftest;
result.ws = load(['/media/transfer/workspace-' str_s 'mod.mat']);
try
[status,cmdout] = system(['ip -f inet address show ' nodeSideNetworkInterface]);
strloc=strfind(cmdout,'inet');
result.computer=cmdout(strloc+5:strloc+16);
catch

end

save([str '_tracking-' num2str(file*offset) '.mat'],'result')
cd(KeyFolder)
fileN = [recordingFolders(scheduleArray(fol,1)).folder '/' recordingFolders(scheduleArray(fol,1)).name str '_tracking-' num2str(file*offset) '.mat'];
save(fileN,'result')



end
cd ..