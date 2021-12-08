%% 
files  = dir(fullfile(input_dir,'*.rhd'));
file_names = {files.name};
file_prefix = ' ';  % optional file prefix before the Intan timestamps 

% Order files for concatentation based on timestamp
% File name format example: 2021-09-06_210906_21062
[dt,di]=sort(datetime(cellfun(@(x) ['20' x(numel(file_prefix)+1:end-4)],file_names,'UniformOutput',false),'InputFormat','yyyy-MM-dd_yyMMdd_HHmmss'));
file_names = file_names(di);

[~,~,frequency_parameters,~] = load_Intan_data_wrapper(files(1).name,aggPath,pwd,DIR(i).name,strsh,selection);
Fs= frequency_parameters.amplifier_sample_rate;
%%
for pre = 1:numel(downhr)
 cd(downloadFolders{downhr(pre)})
DIR=dir('*.rhd');
i=1
[~,~,frequency_parameters,~] = load_Intan_data_wrapper(DIR(i).name,aggPath,pwd,DIR(i).name,strsh,selection);
Fs= frequency_parameters.amplifier_sample_rate;

parfor i= 1:numel(DIR)  % use parfor to load Intan files in session folder in parallel 
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

