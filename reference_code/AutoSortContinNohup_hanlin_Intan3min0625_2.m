% add isSingle, width ,autocorr , calculation/plotting code 
addpath(genpath('/home/luanlab/Documents/MATLAB'))
str_s='06142021'
workspacename= ['workspace-' char(str_s) '.mat']
load(workspacename);
additionalFileName = 'data_' % file namens before the Intan timestamps like 210421_164553
selport='C'
selection = [] % add a list of accepted channels, start from 1, if empty, use all opened channels
mkdir('/media/transfer/tmp')

white =1; 
Folderstr='CMR'

detect_sign=0;
adj_rad=100;
chunk_size = 1E6; % int64;
clip_size=detect_interval*2;
intConvertFactor = 10;% the first is for standard calculation, second is applied on top of converter during  saving data as integer note that the 1/scale 
% is larger than 5 
padding=clip_size*10;
chunk_size_with_padding=chunk_size+2*padding;
aggPath=[aggPath '/Output_RecordControl']
[status,cmdout] = system(['ip -f inet address show ' nodeSideNetworkInterface]);
strloc=strfind(cmdout,'inet');
ipaddress = strip(cmdout(strloc+5:strloc+19),'/')
% to avoid this 
cd /media/transfer
save(['workspace-' str_s 'mod.mat']);



if ChMapNum ==6
load('Mirro_Oversampling_hippo_map.mat');
elseif ChMapNum ==5
load('Oversampling_hippo_map.mat');    
    elseif ChMapNum==7
    load('oversampling_palvo_flex_intanAdapterMap.mat');
    elseif ChMapNum==8
    load('oversampling_palvo_rigid_intanMap.mat');
end

Maps = cell(1,4);
for i =1:4
Maps{i}=Ch_Map_new(:,(1:2)+(i-1)*2);
end
numOfMaps=numel(Maps);



strsh=cell(1,1);
strsh{1}='#!/bin/bash';
strsh{2}='lftp -d -c ''open -e "set ftps:initial-prot ""; \\';
strsh{3}='set ftp:ssl-force true; \\';
strsh{4}='set ftp:ssl-protect-data true; \\';
strsh{5}='open -u hz4888@eid.utexas.edu ftps://ftp.box.com:990; \\';
strsh{6}=['ls ' aggPath ';'''];

fid = fopen('hanlinRun.sh','wt');
for s=1:numel(strsh)
fprintf(fid,strsh{s});
fprintf(fid,'\n');
end
fclose(fid);
%% group files to folders of 30min
cmdout  = dir('*.rhd')
cmdout = {cmdout.name}
r3=cmdout
[dt,di]=sort(datetime(cellfun(@(x) ['20' x(numel(additionalFileName)+1:end-4)],r3,'UniformOutput',false),'InputFormat','yyyyMMdd_HHmmss'));
r3 = r3(di);
currentd = dt(1)
half_an_Hour_Index = zeros(numel(r3),1) % find nearby 31 mins 
init = 1;
while currentd < dt(end)    
    half_an_Hour_Index(find((dt>=currentd) & (dt<currentd + minutes(30))))= init;
    currentd = currentd + minutes(30)
    init = init +1;
end
%% Configure Download 
uniqueHH = unique(half_an_Hour_Index) 
folderCell = arrayfun(@(x) find(half_an_Hour_Index==x), uniqueHH,'UniformOutput',0)
selFolList = cellfun(@(x) x(1),folderCell) % create psudo folder. pick the name of the first file in any 31 mins
downloadFolders = r3(selFolList) % create psudo folder
%% Bulk Transfer
for i = 1:numel(downloadFolders)
mkdir(downloadFolders{i}(1:end-4))
fcCell = folderCell{i};
for fc = 1:numel(fcCell)
movefile(r3{fcCell(fc)},downloadFolders{i}(1:end-4))
end
end
%% decide to download 
if downloadFromBox==1
%% get listing 
[~,cmdout] = unix('bash hanlinRun.sh');
save('cmd','cmdout')
%% analyse listing 
patr=[additionalFileName '\d{6}_\d{6}.rhd'];
r3=regexp(cmdout,patr,'match');
else
cmdout  = dir('**/*.rhd');
cmdout = {cmdout.name};
r3=cmdout
end
[dt,di]=sort(datetime(cellfun(@(x) [yr x(numel(additionalFileName)+1:end-4)],r3,'UniformOutput',false),'InputFormat','yyyyMMdd_HHmmss'));
r3 = r3(di);
currentd = dt(1)
half_an_Hour_Index = zeros(numel(r3),1) % find nearby 31 mins 
init = 1;
while currentd <= dt(end)    

    half_an_Hour_Index(find((dt>=currentd) & (dt<currentd + minutes(30))))= init;
    currentd = currentd + minutes(30)
    init = init +1;
end
%% Configure Download 
numFolderPerBatch = 2;
uniqueHH = unique(half_an_Hour_Index) 
folderCell = arrayfun(@(x) find(half_an_Hour_Index==x), uniqueHH,'UniformOutput',0)
selFolList = cellfun(@(x) x(1),folderCell) % create psudo folder. pick the name of the first file in any 31 mins
pause(wait)
downloadFolders = r3(selFolList); % create psudo folder
downloadFolders = cellfun(@(x) x(1:end-4),downloadFolders,'UniformOutput',false)
numBatch = ceil(numel(downloadFolders)/numFolderPerBatch);
if downloadFromBox==1
else
    for i = 1:numel(downloadFolders)
mkdir(downloadFolders{i})
fcCell = folderCell{i};
% for fc = 1:numel(fcCell)
% movefile(r3{fcCell(fc)},downloadFolders{i})
% end
    end

end
%% Batched Download
 try
parpool('local',1)
 catch
      delete(gcp)
parpool('local',1)
 end

for batch = batchStart:min(numBatch,batchEnd) %--------------------------------------------------------
downfol =  downloadFolders((batch-1)*numFolderPerBatch+1:min(batch*numFolderPerBatch,numel(downloadFolders)))
downhr   = uniqueHH((batch-1)*numFolderPerBatch+1:min(batch*numFolderPerBatch,numel(downloadFolders)))

for pre = 1:numel(downhr)
try 
    cd(downfol{1}) % extremely important switch to check between download from box versus run from local while files are alreayd donwloaded 
    cd ..
catch
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for intan some prepartions

if schedule==1 && batch>1
    downfol=[downloadFolders((batch-1)*numFolderPerBatch) downfol]; % this ensure you start second batch with first batch last folder 
    % and second batch first folder, this happens because you are using 30
    % min segmenations for raw data files but 1hr sorting segments, 30 mins
    % overlapping segments. 
end
%% hack into auto download with below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recordingFolders = dir; % check downloaded folder s 
recordingFolders=recordingFolders([recordingFolders.isdir]);
names = {recordingFolders.name};
nameLength = cellfun(@numel,names);
recordingFolders = recordingFolders(nameLength==(13+numel(additionalFileName))) % check folder sting length 
names = {recordingFolders.name};
[c,ia,ib]=intersect(names,downfol); % check if folder is in the download list 
recordingFolders=recordingFolders(ia);
names=names(ia);Maps



scheduleArray=[];
if schedule==1
scheduleArray(:,1)=1:numel(recordingFolders)-1;
scheduleArray(:,2)=2:numel(recordingFolders);
end
if schedule==2
scheduleArray(:,1)=1:2:numel(recordingFolders);
try
scheduleArray(:,2)=2:2:numel(recordingFolders);    
catch
end
end
scheduleArray

%prompt='Offset by Hr?: ';
offset=1;%input(prompt);
timePerhour = 1; % in hour.


dateStr=str_s;
finalFolder = numOfMaps;mouse ='';pathToDataFile = '';pathToResultFile = '';

totalSession=1:numOfMaps;
files = totalSession;
%% read rhd and convert .mat every 3 mins 

try
parpool('local',4)
catch
delete(gcp)
parpool('local',4)
end

for pre = 1:numel(downhr)
 cd(downloadFolders{downhr(pre)})
DIR=dir('*.rhd');
i=1
[~,~,frequency_parameters,~] = load_Intan_data_wrapper(DIR(i).name,aggPath,pwd,DIR(i).name,strsh,selection);
Fs= frequency_parameters.amplifier_sample_rate;

parfor i= 1:numel(DIR) %parfor
    data=[]
    amplifier_data=[]
    amplifier_channels=[]
    frequency_parameters=[]
    DataSeg = [] 
[amplifier_data,amplifier_channels,frequency_parameters,ADC] = load_Intan_data_wrapper(DIR(i).name,aggPath,pwd,DIR(i).name,strsh,selection);
ports  = {amplifier_channels.port_prefix}
chs_in_this_port  = cell2mat(ports)==selport;
if ~isempty(selection)
chs_in_this_port  = chs_in_this_port(selection)
amplifier_channels=amplifier_channels(selection);
end
amplifier_data=amplifier_data(chs_in_this_port,:);


ADC = int16(ADC*intConvertFactor);
fileN = ['/media/transfer/' DIR(i).name 'ADC.mat'];
parsaveADC(fileN, ADC);


        amplifier_data = single(amplifier_data);


for fol = 1:numOfMaps
mkdir(num2str(fol))
x=[]
x(:,2)=intersect([amplifier_channels.native_order]'+1,Maps{fol});
x(:,1)=1:size(x,1);
Dat_V_Map=x;
cd(num2str(fol))
mkdir(Folderstr)
cd(Folderstr)
parsaveDat_V_Map('Dat_V_Map',Dat_V_Map)
cd ..
cd ..
end

for fol = 1:numOfMaps
    [i fol]
    [c,ia,ib]=intersect([amplifier_channels.native_order]'+1,Maps{fol})
DataSeg = int16(amplifier_data(ia,:)*intConvertFactor);
parsaveDataSeg([num2str(fol) '/' DIR(i).name(1:end-4)],DataSeg)
end


end % finish looping for each file 




unix('rm -rf *.rhd*')
cd .. 
end
delete(gcp)
%% load overlapped hours/or not 

for fol=1:size(scheduleArray,1) % outer loop to go through each segment  folder in this batch
cd([recordingFolders(scheduleArray(fol,1)).name])
KeyFolder = pwd ;
clear x
%
DIR1=dir([recordingFolders(scheduleArray(fol,1)).folder '/' recordingFolders(scheduleArray(fol,1)).name '/*/*.*mat'])
hasNextFol=1;
try
DIR2=dir([recordingFolders(scheduleArray(fol,1)).folder '/' recordingFolders(scheduleArray(fol,2)).name '/*/*.*mat']);
catch
    hasNextFol=0;
end
window_length = Fs;
[b,a]=butter(4,400/(Fs/2),'high');


tic
parpool('local',6)
for sfol = 1:numOfMaps
    
   currentMap = Maps{sfol};
   currentMap = sort(currentMap(:));
%    currentMap =intersect(currentMap,selection);
   Dat_V_Map= importdata([DIR1(1).folder(1:end-2) '/' num2str(sfol) '/' Folderstr '/' 'Dat_V_Map.mat']);
   currentMap =intersect(currentMap,Dat_V_Map(:,2));
   
   if hasNextFol==1
   dataSize  = sum([DIR1.bytes]) +  sum([DIR2.bytes]);
   else
        dataSize  = sum([DIR1.bytes]);
   end
   estimatedLength = dataSize/216e6/128*3600*30E3;
   
   amplifier_data=zeros(numel(currentMap),ceil(estimatedLength),'single');
   
   if hasNextFol==1
   files=[DIR1;DIR2];
   else
       files=DIR1;
   end
   [C,ia,ic]=unique({files.name});
   files = files(ia);
   ep=0;
   
for i=1:numel(files)
    i
    dataS=double(importdata([files(i).folder(1:end-2) '/' num2str(sfol) '/' files(i).name(1:end-4)  '.mat']))';
    dataS = filtfilt(b,a,dataS)/intConvertFactor; % to perform the filter, this must be column vector .
    amplifier_data(:,ep+1:ep+size(dataS,1))= single(dataS)';
    ep=ep+size(dataS,1);
    dataS =[] 
end

dend2=floor(ep/window_length)*window_length;

amplifier_data(:,dend2+1:end)=[];

stdNoise_original = zeros(1,numel(currentMap),'single');

parfor chlfp = 1:  numel(currentMap)
totaldata = single(amplifier_data(chlfp,:));
stdNoise_original(chlfp) = std(totaldata);
end
totaldata=[];
   cd(num2str(sfol))
   cd(Folderstr)
   writemda(int16(intConvertFactor*amplifier_data),['noCMR.mda'],'int16');

  
   
   CMRdotMDA = median(amplifier_data);
   if CMR==1
amplifier_data = bsxfun(@minus,amplifier_data, CMRdotMDA );
elseif CMR==2
    list1=1:15;
    list2=16:31;
amplifier_data(list1,:) = amplifier_data(list1,:) - repmat(median(amplifier_data(list1,:)),size(amplifier_data(list1,:),1),1);
amplifier_data(list2,:) = amplifier_data(list2,:) - repmat(median(amplifier_data(list2,:)),size(amplifier_data(list2,:),1),1);
end
  writemda(int16(intConvertFactor*CMRdotMDA),['CMR.mda'],'int16');
  writemda(int16(amplifier_data(:,1:100000)),[str_s '.mda'],'int16'); %2020 just to make process run 

%% potential whitening 

if white
amplifier_data=amplifier_data';
%     for chstd = 1:size(amplifier_data,1)
%      amplifier_data(chstd,:)= amplifier_data(chstd,:)-mean(amplifier_data(chstd,:));
% end
amplifier_data = bsxfun(@minus,amplifier_data,mean(amplifier_data));

[E,D] = eig(cov(amplifier_data));
amplifier_data=amplifier_data';
[U,S,~]      = svd(amplifier_data,'econ');
amplifier_data=U*(S^-1)*U'*amplifier_data;
amplifier_data=amplifier_data*sqrt(size(amplifier_data,2));
end

   amplifier_data=amplifier_data';
   M = size(amplifier_data,2);
   numCh = M; 
   num_channels=M;
   stdNoise = zeros(1,M);

   t_start = 1;
   t_end = size(amplifier_data,1);
  

   num_timepoints = t_end-t_start+1;
   num_chunks=ceil(num_timepoints/chunk_size);
   N=num_timepoints; 

unix(['rm -rf timeseries.hdf5'])
h5create(['timeseries.hdf5'],'/chunk_size',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/chunk_size', int64(chunk_size))
h5create(['timeseries.hdf5'],'/num_channels',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_channels', int64(num_channels))
h5create(['timeseries.hdf5'],'/num_timepoints',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_timepoints', int64(num_timepoints))
h5create(['timeseries.hdf5'],'/padding',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/padding', int64(padding))
h5create(['timeseries.hdf5'],'/num_chunks',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_chunks', int64(num_chunks))
for m =1:M
           chsig=single(amplifier_data(t_start:t_end,m));
           chsigstd  = chsig;
           pr=prctile(chsig,[25 75]);
           pr25=pr(1);pr75=pr(2);chiqr= pr75-pr25;   
           UL = pr75+1.5*chiqr;LL = pr25-1.5*chiqr;
           chsigstd(chsig>UL)=nan;chsigstd(chsig<LL)=nan; % outlier reject
           stdNoise(m)=nanstd(chsigstd);
           

       for j = 1:num_chunks
           
            padded_chunk=zeros(1,chunk_size_with_padding,'int16') ;
            t1=int64((j-1)*chunk_size); %# first timepoint of the chunk
            t2=int64(min(N,(t1+chunk_size))); %# last timepoint of chunk (+1)
            s1=int64(max(0,t1-padding)); %# first timepoint including the padding
            s2=int64(min(N,t2+padding)); %# last timepoint (+1) including the padding
            aa = padding-(t1-s1)+1;
            if white==1
            padded_chunk(1,aa:aa+s2-s1-1)=int16(chsig((s1+1):s2)*intConvertFactor*intConvertFactor); %# Read the padded chunk
            else
            padded_chunk(1,aa:aa+s2-s1-1)=int16(chsig((s1+1):s2)*intConvertFactor); %# Read the padded chunk
            end
            %[m j sum(isnan(padded_chunk))]
            partName = ['/part-' num2str(m-1) '-' num2str(j-1)];
            h5create(['timeseries.hdf5'],partName,chunk_size_with_padding,'Datatype','int16')
            h5write(['timeseries.hdf5'],partName,padded_chunk)
       end
end
chsig=[];

%%

mode=0;
createGeometryCSVHO;
Dis_mat=pdist2(CSVList,CSVList);
%% perform detection 
save('param','CMR','ChMapNum','stdNoise','stdNoise_original','files') % this variable is later saved to result.param
unix('rm -rf Detect.hdf5')
   for ch = 1:numCh
        data = amplifier_data(t_start:t_end,ch);
        times =  ms4_detect_on_channel(data', detect_threshold*stdNoise(ch),detect_interval,detect_sign,0);
          if numel(times)==0
        times=detect_interval*2+1;
          end
          neighborhood = Dis_mat(ch,:)<=adj_rad;
        val_this_ch = data(times);
     if detect_sign<0
        nearby_neighborhood_maximum0 = movmax(max(-amplifier_data(t_start:t_end,neighborhood),[],2),[detect_interval detect_interval]);
        assign_to_this_neighborhood0=((-val_this_ch)==nearby_neighborhood_maximum0(times));
     elseif detect_sign==0
        nearby_neighborhood_maximum0 = movmax(max(abs(amplifier_data(t_start:t_end,neighborhood)),[],2),[detect_interval detect_interval]);
        assign_to_this_neighborhood0=(abs(val_this_ch)==nearby_neighborhood_maximum0(times));
     elseif detect_sign>0
        nearby_neighborhood_maximum0 = movmax(max(amplifier_data(t_start:t_end,neighborhood),[],2),[detect_interval detect_interval]);
        assign_to_this_neighborhood0=(val_this_ch==nearby_neighborhood_maximum0(times));
     end
      
    validTimes = (times>clip_size/2) & (times < N-clip_size/2);
    times = times(validTimes);
    assign_to_this_neighborhood0 = assign_to_this_neighborhood0(validTimes);
    partName = ['/ch-' num2str(ch) '-times' ];
    h5create('Detect.hdf5',partName,numel(times),'Datatype','int64')
    h5write('Detect.hdf5',partName,int64(times-1))
    partName = ['/ch-' num2str(ch) '-assign' ];
    h5create('Detect.hdf5',partName,numel(assign_to_this_neighborhood0),'Datatype','int8')
    h5write('Detect.hdf5',partName,int8(assign_to_this_neighborhood0))
    sprintf('ch %d detected %d out of %d events ',[ch,sum(assign_to_this_neighborhood0),numel(assign_to_this_neighborhood0)]) 
   end


amplifier_data=[]
savePath=pwd;


cd ..
delete('*.mat')
cd ..
end

%% mountainSort 
try
    delete(gcp)
catch
end
for file=1:numOfMaps % go through each shank
session = [num2str(totalSession(file)) '/CMR']
cd(session)
% unix('mv raw.mda CMR.mda') 
unix(['mv ' str_s '.mda raw.mda']) 
% unix('mv firings.mda firingsWhite.mda') 
mountainSortBash2019
unix(['cp /home/luanlab/Documents/MATLAB/ML_code/ML_code/params.json ./params.json'])
unix(['time ' pythonPath ' sftest.py'])
cd(KeyFolder)
end
files = totalSession;
for file=1:numel(files)
    session = [num2str(totalSession(file)) '/CMR']
cd(session)
param=load('param.mat');
weighted_center_for_mountainSort_large_128EBL('noCMR',session,ChMapNum,intConvertFactor)

% dataFile=str;
  cd(KeyFolder)
end


%% Configure Upload
     str = 'noCMR';

     
     
for file=1:numel(files)
    session = [num2str(file) '/' Folderstr];

     cd([pathToDataFile mouse session ])
try
    unix('rm -f noCMR.mda')
    unix('rm -f timeseries.hdf5')
    unix('rm -f CMR.mda')
catch
end

msp = loadjson(['params.json']);
sftest = fileread('sftest.py');



folderStru = dir(['**/' str '_tracking.mat']);
Sorting_result = load([folderStru.folder '/' str '_tracking.mat']);

Sorting_result =Sorting_result.result;

for unit=1:numel(Sorting_result.list)
wav = Sorting_result.waveform{unit};
wav_std = Sorting_result.waveformStd{unit};

P2P = min(wav');
ch=find(P2P==min(P2P));
ch=ch(1);
tp=find(wav(ch,:)==min(P2P));
tp=tp(1);
Sorting_result.peakCV(unit)=abs(wav_std(ch,tp)/min(P2P));

FiringTimeForThisUnit = Sorting_result.time{unit};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/(Fs/1E3); % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
maxTime = max(cell2mat(Sorting_result.time))/Fs;

Sec_bins = double(FiringTimeForThisUnit)/Fs;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;


Sorting_result.Max_Ins_5s_FR(unit)=Max_Ins_5s_FR;
Sorting_result.Avg_FR(unit)=numel(FiringTimeForThisUnit)/maxTime;

try
Sorting_result.multi1(unit)= counts(1);
Sorting_result.multi2(unit)= counts(2);
Sorting_result.single(unit)= ((counts(2)+counts(1))/numel(FiringTimeForThisUnit))<=0.01;
Sorting_result.violation(unit)= ((counts(2)+counts(1))/numel(FiringTimeForThisUnit));
catch
    try
Sorting_result.multi2(unit)= counts(2);
    catch
    end
end
% we allow one spike per min.
end
result = Sorting_result;
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
end
end

unix('rm -rf /tmp/*MountainSort*')
unix('rm -rf /media/transfer/tmp/*')
%% plotting 
parpool('local',16)
layernum=[1 2 3 4];
HD=0;
CodeForPlottingOnLocalWindows32Channel
