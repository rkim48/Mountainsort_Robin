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

selection = []; % list of accepted channels (start from index 1), if empty, use all channels
load('pavlo_1x32_intanMap_v2.mat'); % 1x32 linear design
%%

concat_data = [];
date_arr = '';
parfor i = 1:numel(file_names)  % use parfor to load Intan files in session folder in parallel 
    [amplifier_data,amplifier_channels,frequency_parameters,~] = load_Intan_data_wrapper(fullfile(input_dir,file_names{i}),selection);
    concat_data = [concat_data, int16(amplifier_data*10)]; % multiply by 10 to save 1st decimal place when converting to int
end

Fs= frequency_parameters.amplifier_sample_rate;

writemda(concat_data,fullfile(input_dir,'data.mda'),'int16');

cd(output_dir)

%% Change geom.csv based on recorded channels

[~,~,good_ch] = intersect([amplifier_channels.native_order]',sort(channel_map));
map = generateMap(channel_map,good_ch,input_dir);

