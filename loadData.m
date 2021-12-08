%% 
addpath(input_dir)
addpath('parsave_utils')
files  = dir(fullfile(input_dir,'*.rhd'));
file_names = {files.name};
file_prefix = ' ';  % optional file prefix before the Intan timestamps 

% Order files for concatentation based on timestamp
% File name format example: 2021-09-06_210906_21062
[dt,di]=sort(datetime(cellfun(@(x) ['20' x(numel(file_prefix)+1:end-4)],file_names,'UniformOutput',false),'InputFormat','yyyy-MM-dd_yyMMdd_HHmmss'));
file_names = file_names(di);

[~,~,frequency_parameters,~] = load_Intan_data_wrapper(files(1).name,aggPath,pwd,DIR(i).name,strsh,selection);
Fs= frequency_parameters.amplifier_sample_rate;
selection = [];
load('pavlo_1x32_intanMap_v2.mat'); % 1x32 linear design
%%

concat_data = [];
date_arr = '';
parfor i = 1:numel(file_names)  % use parfor to load Intan files in session folder in parallel 
    [amplifier_data,amplifier_channels,frequency_parameters,~] = load_Intan_data_wrapper(fullfile(input_dir,file_names{i}),selection);

    
    x2 =intersect([amplifier_channels.native_order]',channel_map);
    x1=1:numel(x2);
    Dat_V_Map=[x1' x2];
    parsaveDat_V_Map('Dat_V_Map',Dat_V_Map) % TO-DO: save in specific directory
    
    
    concat_data = [concat_data, int16(amplifier_data*10)]; % multiply by 10 to save 1st decimal place when converting to int
%     parsaveDataSeg([output_dir '/' file_names{i}(1:end-4)],single(amplifier_data)); % convert to single precision  
    
end
Fs= frequency_parameters.amplifier_sample_rate;

writemda(concat_data,fullfile(input_dir,'data.mda'),'int16');

cd(output_dir)

