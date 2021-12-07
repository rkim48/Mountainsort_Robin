% addpath(genpath('/home/luanlab/Documents/MATLAB'))
workspacename = 'workspace.mat';
load(workspacename);

selport = 'C'; % port of interest 
selection = [] % add a list of accepted channels, start from 1, if empty, use all opened channels
mkdir('temp')
white =1; 
Folderstr = 'CMR'
detect_sign=0;
adj_rad=100;
chunk_size = 1E6; % int64;
clip_size=detect_interval*2;
intConvertFactor = 10; % the first is for standard calculation, second is applied on top of converter during saving data as integer note that the 1/scale 
% is larger than 5 
padding=clip_size*10;
chunk_size_with_padding=chunk_size+2*padding;
ChMapNum = 2;
good_channels = []; % add a list of accepted channels, start from 1, if empty, use all opened channels
pythonPath = '/home/robin/miniconda3/envs/mountainlab/bin/python';


source_dir = uigetdir;
save_dir = uigetdir;
mkdir transfer

%% Channel map selection
if ChMapNum ==1
    load('oversampling_palvo_flex_intanAdapterMap.mat');
    Maps = cell(1,4);
    for i =1:4
        Maps{i}=Ch_Map_new(:,(1:2)+(i-1)*2);
    end
    nMAPS=numel(Maps);
elseif ChMapNum==2
    load('pavlo_1x32_intanMap_v2.mat');
    Maps = cell(1,1);
    nMAPS = 1;
    Maps{1} = Ch_Map_new;
end



%% We want to read files locally 

rhd_files = dir(fullfile(source_dir,'*.rhd')); % list all files with .rhd extension

rhd_name = rhd_files(1).name;
rhd_path = fullfile(rhd_files(1).folder,rhd_name);
[amplifier_data,amplifier_channels,frequency_parameters,ADC]=read_Intan_RHD2000_fileV3(rhd_path,good_channels);
Fs = frequency_parameters.amplifier_sample_rate;
    
parfor i = 1:numel(rhd_files)
    rhd_name = rhd_files(i).name;
    rhd_path = fullfile(rhd_files(i).folder,rhd_name);
    [amplifier_data,amplifier_channels,frequency_parameters,ADC]=read_Intan_RHD2000_fileV3(rhd_path,good_channels);
    ports  = {amplifier_channels.port_prefix};
    chs_in_this_port  = cell2mat(ports) == selport;
    if ~isempty(good_channels)
        chs_in_this_port  = chs_in_this_port(good_channels)
        amplifier_channels=amplifier_channels(good_channels);
    end
    amplifier_data=amplifier_data(chs_in_this_port,:);
    ADC = int16(ADC*intConvertFactor);
    
    fileN = ['transfer/' rhd_name 'ADC.mat']; % why save adc.mat?
    parsaveADC(fileN, ADC);
    
    amplifier_data = single(amplifier_data);
    
    for folder = 1:nMAPS
        mkdir(num2str(folder))
        x =intersect([amplifier_channels.native_order]'+1,Maps{folder});
        x =[(1:numel(x))' x];
        Dat_V_Map=x;
        cd(num2str(folder))
        mkdir(Folderstr)
        cd(Folderstr)
        parsaveDat_V_Map('Dat_V_Map',Dat_V_Map)
        cd ../..
    end
    
    for folder = 1:nMAPS
        [c,ia,ib]=intersect([amplifier_channels.native_order]'+1,Maps{folder})
        data_seg = int16(amplifier_data(ia,:)*intConvertFactor);
        parsaveDataSeg([num2str(folder) '/' extractBefore(rhd_name,'.rhd')],data_seg)
    end

end

%%

window_length = Fs;
[b,a]=butter(4,400/(Fs/2),'high');

files=dir(source_dir);
files = files(~cellfun(@(c)isequal(c,0), {files.bytes}));
for folder = 1:nMAPS
   currentMap = Maps{sfol};
   currentMap = sort(currentMap(:));
   
   Dat_V_Map= importdata([num2str(folder) '/' Folderstr '/' 'Dat_V_Map.mat']);
   currentMap =intersect(currentMap,Dat_V_Map(:,2));
   
   dataSize  = sum([files.bytes]);

   estimatedLength = dataSize/216e6/128*3600*30E3;
   
   amplifier_data=zeros(numel(currentMap),ceil(estimatedLength),'single');
   
   [C,ia,ic]=unique({files.name});
   files = files(ia);
   ep=0;
   
for i=1:numel(files)
    dataS=double(importdata([num2str(folder) '/' files(i).name(1:end-4) '.mat']))';
    dataS = filtfilt(b,a,dataS)/intConvertFactor; % to perform the filter, this must be column vector .
    amplifier_data(:,ep+1:ep+size(dataS,1))= single(dataS)';
    ep=ep+size(dataS,1);
    dataS = [] ;
end

dend2=floor(ep/window_length)*window_length;
amplifier_data(:,dend2+1:end)=[];
stdNoise_original = zeros(1,numel(currentMap),'single');
parfor chlfp = 1: numel(currentMap)
    totaldata = single(amplifier_data(chlfp,:));
    stdNoise_original(chlfp) = std(totaldata);
end

totaldata=[];
cd(num2str(folder))
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
writemda(int16(amplifier_data(:,1:100000)),[str_s '.mda'],'int16');  %2020 just to make process run 

if white
    amplifier_data=amplifier_data';
    amplifier_data = bsxfun(@minus,amplifier_data,mean(amplifier_data));
    [E,D] = eig(cov(amplifier_data));
    amplifier_data=amplifier_data';
    [U,S,~]  = svd(amplifier_data,'econ');
    amplifier_data=U*(S^-1)*U'*amplifier_data;
    amplifier_data=amplifier_data*sqrt(size(amplifier_data,2));
end

amplifier_data=amplifier_data';
nCH = size(amplifier_data,2);
stdNoise = zeros(1,nCH);

t_start = 1;
t_end = size(amplifier_data,1);
nCHUNKS=ceil(num_timepoints/chunk_size);
nTIMEPTS =t_end-t_start+1;

unix(['rm -rf timeseries.hdf5'])
h5create(['timeseries.hdf5'],'/chunk_size',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/chunk_size', int64(chunk_size))
h5create(['timeseries.hdf5'],'/num_channels',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_channels', int64(nCH))
h5create(['timeseries.hdf5'],'/num_timepoints',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_timepoints', int64(nTIMEPTS))
h5create(['timeseries.hdf5'],'/padding',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/padding', int64(padding))
h5create(['timeseries.hdf5'],'/num_chunks',1,'Datatype','int64')
h5write(['timeseries.hdf5'], '/num_chunks', int64(nCHUNKS))

for m =1:nCH
    ch_sig=single(amplifier_data(t_start:t_end,m));
    ch_sig_std  = ch_sig;
    pr=prctile(ch_sig,[25 75]);
    pr25=pr(1);pr75=pr(2);chiqr= pr75-pr25;   
    UL = pr75+1.5*chiqr;LL = pr25-1.5*chiqr;
    ch_sig_std(ch_sig>UL)=nan;ch_sig_std(ch_sig<LL)=nan; % outlier reject
    stdNoise(m)=nanstd(ch_sig_std);
           
   for j = 1:num_chunks
        padded_chunk=zeros(1,chunk_size_with_padding,'int16') ;
        t1=int64((j-1)*chunk_size); %# first timepoint of the chunk
        t2=int64(min(N,(t1+chunk_size))); %# last timepoint of chunk (+1)
        s1=int64(max(0,t1-padding)); %# first timepoint including the padding
        s2=int64(min(N,t2+padding)); %# last timepoint (+1) including the padding
        aa = padding-(t1-s1)+1;
            if white==1
                padded_chunk(1,aa:aa+s2-s1-1)=int16(ch_sig((s1+1):s2)*intConvertFactor*intConvertFactor); %# Read the padded chunk
            else
                padded_chunk(1,aa:aa+s2-s1-1)=int16(ch_sig((s1+1):s2)*intConvertFactor); %# Read the padded chunk
            end
        partName = ['/part-' num2str(m-1) '-' num2str(j-1)];
        h5create(['timeseries.hdf5'],partName,chunk_size_with_padding,'Datatype','int16')
        h5write(['timeseries.hdf5'],partName,padded_chunk)
   end
end

ch_sig=[];
mode=0;
createGeometryCSVHO2;
Dis_mat=pdist2(CSVList,CSVList);

save('param','CMR','ChMapNum','stdNoise','stdNoise_original','files') % this variable is later saved to result.param
unix('rm -rf Detect.hdf5')

for ch = 1:nCH
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

for file=1:nMAPS% go through each shank
    session = [num2str(file) '/CMR']
    cd(session)
    % unix('mv raw.mda CMR.mda') 
    unix(['mv ' str_s '.mda raw.mda']) 
    % unix('mv firings.mda firingsWhite.mda') 
    
    
    mountainSortBash2019 % Needs to be fixed
    
    
    
    
    unix(['cp /media/robin/Shared/ML_code/params.json ./params.json'])
    unix(['time ' pythonPath ' sftest.py'])
    cd(KeyFolder)
end

files = totalSession;
for file=1:numel(files)
    session = [num2str(file) '/CMR']
    cd(session)
    param=load('param.mat');
    weighted_center_for_mountainSort_large_128EBL('noCMR',session,ChMapNum,intConvertFactor)
    cd(KeyFolder)
end

