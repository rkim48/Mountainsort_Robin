% Setup file
%% Use conda m-file to switch to conda environment with Mountainsort code
% Initialize
conda.init
% List available conda environments
conda.getenv
% Set environment 
% e.g. conda.setenv('mountainlab2')
conda.setenv('mountainlab2')
%% Change ML config variables if needed (ignore if no error) 
% Do not set output directory to server due to hard link issues. Set output
% directory to local machine

% system('ml-config')
% system('ML_PACKAGE_SEARCH_DIRECTORY=/home/robin/miniconda3/envs/mountainlab2/etc/mountainlab/packages')
% system('ML_CONFIG_FILE=/home/robin/miniconda3/envs/mountainlab2/etc/mountainlab/mountainlab.env')
%% Mountainsort parameters
M = MountainsortParams();

working_dir =  M.working_dir;
copy_str = M.copy_str;
input_dir = M.input_dir;
output_dir = M.output_dir;
file_name = M.file_name;
samplerate = M.samplerate;
geom_file = M.geom_file;
clip_size = M.clip_size;
freq_min = M.freq_min;
freq_max = M.freq_max;
detect_interval = M.detect_interval;
detect_threshold = M.detect_threshold;
adjacency_radius = M.adjacency_radius;
detect_sign = M.detect_sign;

% Load channel map and sort Intan files by date
impedance_threshold = 2e6; % accept channels under 2 MOhm impedance

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
% get characters between '_' and '.'
file_timestamps = regexp({files.name}, '(?<=_)[^-]+(?=\.)', 'match', 'once'); 
% sort based on timestamps
[dt,di]=sort(datetime(cellfun(@(x) x,file_timestamps,'UniformOutput',false),'InputFormat','yyMMdd_HHmmss'));
file_names = file_names(di);

% Concatenate data 
concat_data=[];
parfor i = 1:numel(file_names)  % use parfor to load Intan files in session folder in parallel 
    if rhs == 1
        [amplifier_data,amplifier_channels,frequency_parameters,channel_struct] = ...
            load_IntanRHS_data_wrapper(fullfile(input_dir,file_names{i}));
    else
        [amplifier_data,amplifier_channels,frequency_parameters,channel_struct] = ...
            load_IntanRHD_data_wrapper(fullfile(input_dir,file_names{i}));
    end
    parsaveIntanStruct(amplifier_channels,frequency_parameters,working_dir);
    concat_data = [concat_data, amplifier_data]; 
end

load('Intan_struct');

all_impedances = [Intan_struct.amplifier_channels.electrode_impedance_magnitude];
native_order = [Intan_struct.amplifier_channels.native_order];
low_impedance_idx = find(all_impedances < impedance_threshold);
low_impedance_ch = native_order(low_impedance_idx);

Fs= Intan_struct.frequency_parameters.amplifier_sample_rate;
samplerate = Fs;
desired_length = ceil(length(concat_data)/Fs)*Fs;
nZEROS = desired_length - length(concat_data);
concat_data = concat_data(low_impedance_idx,:); 

concat_data = concat_data - mean(concat_data);
concat_data_appended = [concat_data zeros(size(concat_data,1),nZEROS)];
% multiply by 10 to save 1st decimal place when converting to int16
writemda(int16(10*concat_data_appended),fullfile(input_dir,'data.mda'),'int16');

clear concat_data concat_data_appended amplifier_data
% Change geom.csv based on recorded channels
cd(working_dir)
[~,~,good_ch] = intersect(low_impedance_ch,sort(channel_map));
[map, depth_idx] = generateMap(channel_map,good_ch,input_dir);
% Channel index table 
% col 1 - 1:nCH where nCH is number of good channels 
% col 2 - depth index (1 shallow, 32 deep)
% col 3 - native Intan order (0-31)
col1 = 1:numel(good_ch);
col2 = depth_idx;
col3 = low_impedance_ch; 
ch_idx_table = array2table([col1' col2 col3'],'VariableNames',{'1:nCH','Depth index','Native Intan order'});