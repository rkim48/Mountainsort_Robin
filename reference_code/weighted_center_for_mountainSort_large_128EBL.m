function weighted_center_for_mountainSort_large_128EBL(dataFile,path2file,ChMapNum,intConvertFactor)
% ch_num= length(amplifier_channels);
% figure(1);clf; co=get(gca,'colororder');
% offset= max(max(amplifier_data))*0.5;

% d = fdesign.highpass('Fst,Fp,Ast,Ap',0.015,0.025,60,1);
%    Hd = design(d,'equiripple');
%  0    -0?=-01234509??
% %    for i=1:ch_num
% i=19;
%        datafilt(i,:)  = filter(Hd,amplifier_data(i,:));
%        figure;
%        plot(t_amplifier(1:1447200), datafilt(i,1:1447200)+(i-1)*300,'-','color',co(mod(i,7)+1,:)); hold on;
%        text(max(t_amplifier)+2,(i-1)*offset,sprintf(amplifier_channels(i).native_channel_name),'color',co(mod(i,7)+1,:));
% %    end3
%
% ok ,we need the reading thing to get Channel Maps.
% clear all
warning off

load('Dat_V_Map.mat')


m = memmapfile([dataFile '.mda'],...
    'Offset',20,...
    'Format','int16',...
    'Writable',true);
nChansInRawFile=size(Dat_V_Map,1);
numOfCol = numel(m.data)/nChansInRawFile;

str=dataFile;
Source = memmapfile([path2file '/' dataFile '.mda'],'Offset',20, 'Format', {'int16', [nChansInRawFile, numOfCol],'x'},'Writable',true);
%
% end





if ChMapNum ==1
    load('Ebeam32by4Map.mat'); % for non-06252019 mice
elseif ChMapNum==0
    load('ch_map_pink.mat');
elseif ChMapNum==4
    load('D:\Box Sync\0919 stroke\Finger_map.mat');
elseif ChMapNum==3
    load('D:\0925 stroke\recording 2017-10-11\PCB_Finger.mat');
elseif ChMapNum==2
    load('NewMap_Hippo_large.mat');
elseif ChMapNum==5
    load('Oversampling_hippo_map.mat');
elseif ChMapNum==6
    load('Mirro_Oversampling_hippo_map.mat');
    elseif ChMapNum==7
    load('oversampling_palvo_flex_intanAdapterMap.mat');
    elseif ChMapNum==8
    load('oversampling_palvo_rigid_intanMap.mat');
    elseif ChMapNum==9
    load('palvo_1x32_intanMap.mat');
    
end
Ch_Map=Ch_Map_new;
if ChMapNum==5 || ChMapNum==6 || ChMapNum==7 || ChMapNum==8 
    intsecNum=cellfun(@(x) numel(intersect(Dat_V_Map,x)),Maps);
    [a,b]=max(intsecNum);
    Ch_Map=Maps{b};
    Ch_Map=[flip(Ch_Map(1:8,:)');flip(Ch_Map(9:16,:)')]; % for later plotting, I know this is strange opt
    Ch_Map_new=Ch_Map;
end

for i=1:size(Ch_Map,1)
    for j=1:size(Ch_Map,2)
        if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j), 1))
            Ch_Map(i,j)=0;
        end
    end
end

Ch_Map_2 = Ch_Map;
for i=1:size(Ch_Map,1)
    for j=1:size(Ch_Map,2)
        if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
            Ch_Map_2(i,j)=0;
        else
            Ch_Map_2(i,j)=find(Dat_V_Map(:,2)==Ch_Map(i,j));
        end
    end
end




nChansInRawFile=size(Dat_V_Map,1);
nChansInDatFile=size(Dat_V_Map,1);
tBefore=30;
tAfter=45;



ElecList=dir('*.csvn'); % get all electrical recording spike timing data file spec.
ImpedanceRec=zeros(32,size(ElecList,1));



for i=1:size(ElecList,1)
    book=importdata(ElecList(i).name);
    
    impedance=book.data(:,2);
    
    if numel(impedance)>32
        ImpedanceRec(:,i)=impedance(17:48);
    else
        ImpedanceRec(:,i)=impedance;
    end
    
end



if size(ImpedanceRec,2)~=1
    AvgImpedance=median(ImpedanceRec')';
else
    AvgImpedance=ImpedanceRec;
end

%AvgImpedance(AvgImpedance>2E6)=2E6;

Ch_Map_Imp=Ch_Map_new-15;
for row=1:4
    for col=1:8
        newImpMatrix(row, col) = 0;% AvgImpedance(Ch_Map_Imp(row,col));
        if newImpMatrix(row, col)>3E6
            newImpMatrix(row, col)= 3E6;
        end
    end
end
newImpMatrix=newImpMatrix/1E6;
% cd(currentpath)
CMRdotMDA = readmda('CMR.mda',1);

mdaStrut = dir('**/firings.mda');
mdaName =[mdaStrut.folder '/firings.mda'];
% [file, path, ~] = ...
% uigetfile('*.kwik', 'Select an KWIK', 'MultiSelect', 'off');
% KwikName=[path,file];%dir('*.kwik');
A=readmda(mdaName,0);

cluster=A(3,:);
MyTimes=A(2,:)+1; % handle difference between python and matlab
priCh = A(1,:)+1;
cluterPriCh=unique([cluster' priCh'],'rows');



color=jet(120);
color=color(5:104,:);
firstRow=[0 0 0];
color=[firstRow;color];

valid = MyTimes>(tBefore) &  MyTimes<(numOfCol-tAfter);
MyTimes=MyTimes(valid);
cluster=cluster(valid);
% Mask=Mask(valid,:);
% valid= (floor(cluster)==cluster)&(cluster<=numel(unique(cluster)));
% MyTimes=MyTimes(valid);
% cluster=cluster(valid);

list=unique(cluster);


pdfCellArray=cell(1,size(list,1));

WCList = zeros(numel(list),3);

% Prepare saving structure
mkdir('Refine_mountain')
cd ('Refine_mountain')
% oldList = cell2mat(raw(:,1));
% newList = zeros(size(oldList));

result.date=str;
%             result.session = folderNum;
result.Dat_V_Map=Dat_V_Map;
result.Ch_Map=Ch_Map;
result.Ch_Map_2=Ch_Map_2;
result.selectedClusters=list;

result.waveform=cell(size(result.selectedClusters));
result.waveformStd=cell(size(result.selectedClusters));
%             result.Intensity=cell(size(result.selectedClusters));
result.P2P = cell(size(result.selectedClusters));
result.Location = cell(size(result.selectedClusters));
result.ValleyAcrossTime = cell(size(result.selectedClusters));
result.time =cell(size(result.selectedClusters));
result.Max_Ins_5s_FR=zeros(size(result.selectedClusters));
result.Avg_FR=zeros(size(result.selectedClusters));
result.newImpMatrix=newImpMatrix;
result.AvgImpedance=AvgImpedance;
result.PeakAcrossTime = cell(size(result.selectedClusters));
%
SecValley = cell(1,numel(list));

for clu=1:numel(list)%[3 57]
    
    newList(clu)=1;
    
    
    
    clu
    PrimaryChannel = find(cluterPriCh(:,1)==list(clu));
    PrimaryChannel=cluterPriCh(PrimaryChannel,2);
    SelectedPriLoc = zeros(size(Ch_Map,1),size(Ch_Map,2));
    [row col] = find(Ch_Map_2==PrimaryChannel);
    SelectedPriLoc(row,col)=1;
    
    
    pdfCellArray{clu}=[num2str(list(clu)) '.pdf'];
    p=find(cluster==list(clu));
    FiringTimeForThisUnit = MyTimes(p);
    result.time{clu}=MyTimes(p);
    
    ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds.
    [counts] = histc(ISI_in_MS_bins,0:1:50);
    [counts_M] = histc(ISI_in_MS_bins,0:20:1000);
    [counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
    
    
    Sec_bins = double(FiringTimeForThisUnit)/20E3;
    [counts_5s_1] = histc(Sec_bins,0:5:720);
    [counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
    Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
    
    
    
    unit=[]; %MeanPlot
    unitCM=[];
    unit=zeros(nChansInDatFile,tBefore+tAfter+1,length(p),'single');
    unitCM=zeros(nChansInDatFile,tBefore+tAfter+1,length(p),'single');
    
    for k=1:length(p)
        temp=Source.Data.x(1:nChansInDatFile,MyTimes(p(k))-tBefore:MyTimes(p(k))+tAfter);
        unit(:,:,k)=1/intConvertFactor*single(temp);
        unitCM(:,:,k)=1/intConvertFactor*single(bsxfun(@minus,temp,CMRdotMDA(MyTimes(p(k))-tBefore:MyTimes(p(k))+tAfter)));
        
    end
    % ValleyAcrossTime = min(reshape(unit,[size(unit,1)*size(unit,2),size(unit,3)]));
    
    MeanPlot=mean(unit,3);
    if size(unit,3)>3.6E6
        StdPlot=std(unit(:,:,1:1E5),0,3);
    else
        StdPlot=std(unit,0,3);
    end
    P2P = max(MeanPlot')-min(MeanPlot');
    ch=find(P2P==max(P2P));
    ch=ch(1);
    [sortP2P,sP2P]=sort(P2P,'descend');
    PeakAcrossTime = unit(PrimaryChannel+1,tBefore+1,:);
    ValleyAcrossTime = unit(ch,tBefore+1,:);
    ValleyAcrossTime  =ValleyAcrossTime(:);
    PeakAcrossTime    =PeakAcrossTime(:);
    
    result.waveform{clu}=MeanPlot;
    result.P2P{clu}=P2P;
    result.waveformStd{clu}=StdPlot;
    maxV = max(abs(PeakAcrossTime),abs(PeakAcrossTime));
    if maxV<3276
    result.ValleyAcrossTime{clu}=int16(ValleyAcrossTime);
    result.PeakAcrossTime{clu}=int16(PeakAcrossTime);
    else
    result.ValleyAcrossTime{clu}=ValleyAcrossTime;
    result.PeakAcrossTime{clu}=PeakAcrossTime;
    end
    
    result.Location{clu}=WCList(clu,:);
    result.Max_Ins_5s_FR(clu)= Max_Ins_5s_FR;
    result.Avg_FR(clu)= sum(cluster==list(clu))/(numOfCol/20E3);
    %--------------------------------------------repeat 
    unit=unitCM;
    MeanPlot=mean(unit,3);
    if size(unit,3)>3.6E6
        StdPlot=std(unit(:,:,1:1E5),0,3);
    else
        StdPlot=std(unit,0,3);
    end
    P2P = max(MeanPlot')-min(MeanPlot');
    ch=find(P2P==max(P2P));
    ch=ch(1);
    [sortP2P,sP2P]=sort(P2P,'descend');
    PeakAcrossTime = unit(PrimaryChannel,tBefore+1,:);
    ValleyAcrossTime = unit(ch,tBefore+1,:);
    ValleyAcrossTime  =ValleyAcrossTime(:);
    PeakAcrossTime    =PeakAcrossTime(:);
    
    maxV = max(abs(PeakAcrossTime),abs(PeakAcrossTime));
    
    result.waveform_CM{clu}=MeanPlot;
    result.P2P_CM{clu}=P2P;
    result.waveformStd_CM{clu}=StdPlot;
    if maxV<3276
    result.ValleyAcrossTime_CM{clu}=int16(ValleyAcrossTime);
    result.PeakAcrossTime_CM{clu}=int16(PeakAcrossTime);
    else
    result.ValleyAcrossTime_CM{clu}=ValleyAcrossTime;
    result.PeakAcrossTime_CM{clu}=PeakAcrossTime;
    end
end
result.list=list;
result.clusterPriCh=cluterPriCh;

save([result.date '_tracking'],'result')
