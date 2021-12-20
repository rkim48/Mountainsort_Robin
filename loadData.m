%%
working_dir =  '/media/robin/Shared/mountainsort_Robin';
input_dir = uigetdir;
output_dir = fullfile(input_dir,'outputs'); mkdir(output_dir);
samplerate=20000;
geom_file = 'geom.csv';
clip_size=60;
freq_min=300;
freq_max=5000;
detect_interval=27;
detect_threshold=4.5;
adjacency_radius=100;
detect_sign=-1;
%% 
addpath(working_dir)
addpath(input_dir)
files  = dir(fullfile(input_dir,'*.rhs'));

if ~isempty(files) 
    rhs = 1;
    load('pavlo_1x32_intanMap_RHS.mat'); % 1x32 linear design
else
    rhs = 0;
    files  = dir(fullfile(input_dir,'*.rhd'));
    load('pavlo_1x32_intanMap_RHD.mat'); % 1x32 linear design
end

file_names = {files.name};
file_prefix = '';  % optional file prefix before the Intan timestamps 

% Order files for concatentation based on timestamp
% File name format example: 2021-09-06_210906_21062
% (yyyy-MM-dd_yyMMdd_HHmmss)
% OR 211210_165550 (yyMMdd_HHmmss)

[dt,di]=sort(datetime(cellfun(@(x) ['20' x(numel(file_prefix)+1:end-4)],file_names,...
    'UniformOutput',false),'InputFormat','*yyMMdd_HHmmss'));

% [dt,di]=sort(datetime(cellfun(@(x) [x(numel(file_prefix)+1:end-4)],file_names,'UniformOutput',false),...
%     'InputFormat','yyMMdd_HHmmss'));
file_names = file_names(di);

selection = []; % list of accepted channels (start from index 1), if empty, use all channels

%%
concat_data=[];
parfor i = 1:numel(file_names)  % use parfor to load Intan files in session folder in parallel 
    if rhs == 1

        [amplifier_data,amplifier_channels,frequency_parameters,~] = load_IntanRHS_data_wrapper(fullfile(input_dir,file_names{i}),selection);
    else
        
        [amplifier_data,amplifier_channels,frequency_parameters,~] = load_IntanRHD_data_wrapper(fullfile(input_dir,file_names{i}),selection);
    end
    concat_data = [concat_data, int16(amplifier_data*10)]; % multiply by 10 to save 1st decimal place when converting to int
end

Fs= frequency_parameters.amplifier_sample_rate;
desired_length = ceil(length(concat_data)/Fs)*Fs;
nZEROS = desired_length - length(concat_data);
concat_data_appended = [concat_data zeros(size(concat_data,1),nZEROS)];

writemda(concat_data_appended,fullfile(input_dir,'data.mda'),'int16');

cd(output_dir)
clear concat data concat_data_appended amplifier_data

files = dir(fullfile(input_dir,'*.mda')); file_names ={files.name};
file_name = file_names{1};

% Change geom.csv based on recorded channels
cd(working_dir)
[~,~,good_ch] = intersect([amplifier_channels.native_order]',sort(channel_map));
map = generateMap(channel_map,good_ch,input_dir);

