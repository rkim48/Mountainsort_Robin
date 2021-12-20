% Setup file
addpath('mdaio/')
addpath('channel_maps/')
addpath('util/')
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
%% Specify sorting arguments 
working_dir =  '/media/robin/Shared/mountainsort_Robin'; % change this
input_dir = uigetdir;
output_dir = fullfile(input_dir,'outputs'); mkdir(output_dir); % saves outputs to input_dir/outputs
% output_dir = uigetdir;
samplerate=20000;
geom_file = 'geom.csv';
clip_size=60;
freq_min=300;
freq_max=5000;
detect_interval=27;
detect_threshold=4.5;
adjacency_radius=100;
detect_sign=-1;
%% Pass arguments to copy of Mountainsort bash script 
addpath(working_dir)
cd(working_dir)
bash_script = 'Mountainsort_Robin.sh';
execute_dir = working_dir; % where the bash script copy is saved and where it is executed from 
copy_str = strcat(extractBefore(bash_script,'.sh'),'_copy.sh'); 
bash_script_copy = fullfile(execute_dir,copy_str); % make copy of original .sh file

copyfile(bash_script,bash_script_copy);  % copy and move to execute directory 
fileattrib(bash_script_copy,'+w')

% Read in MountainSort .sh file and place each line in a cell
fid = fopen(bash_script_copy,'r');
nn=1;
tline = fgetl(fid);
A{nn} = tline;

while ischar(tline)
    nn = nn+1;
    tline = fgetl(fid);
    A{nn} = tline;
end

fclose(fid);

% Write arguments to copied .sh file
A{1} = ['input_dir=',input_dir];
A{2} = ['output_dir=',output_dir];
A{3} = ['samplerate=',num2str(samplerate)];
A{4} = ['file_name=',file_name];
A{5} = ['geom_file=',geom_file];
A{6} = ['clip_size=',num2str(clip_size)];
A{7} = ['freq_min=',num2str(freq_min)];
A{8} = ['freq_max=',num2str(freq_max)];
A{9} = ['detect_interval=',num2str(detect_interval)];
A{10} = ['detect_threshold=',num2str(detect_threshold)];
A{11} = ['adjacency_radius=',num2str(adjacency_radius)];
A{12} = ['detect_sign=',num2str(detect_sign)];

fid = fopen(bash_script_copy, 'w');
for nn = 1:numel(A)
    if A{nn+1} == -1
        fprintf(fid,'%s', A{nn});
        break
    else
        fprintf(fid,'%s\n', A{nn});
    end
end
fclose(fid);

fileattrib(bash_script_copy,'+x')
%% Load channel map and sort Intan files by date
selection = [];  % list of accepted channels (start from index 1), if empty, use all channels
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
%% Concatenate data 
concat_data=[];
parfor i = 1:numel(file_names)  % use parfor to load Intan files in session folder in parallel 
    if rhs == 1
        [amplifier_data,amplifier_channels,frequency_parameters] = ...
            load_IntanRHS_data_wrapper(fullfile(input_dir,file_names{i}),selection);
    else
        [amplifier_data,amplifier_channels,frequency_parameters] = ...
            load_IntanRHD_data_wrapper(fullfile(input_dir,file_names{i}),selection);
    end
    parsaveIntanStruct(amplifier_channels,frequency_parameters,workicd ..ng_dir);
    % multiply by 10 to save 1st decimal place when converting to int
    concat_data = [concat_data, int16(amplifier_data*10)]; 
end

load('Intan_struct');

all_impedances = [Intan_struct.amplifier_channels.electrode_impedance_magnitude];
native_order = [Intan_struct.amplifier_channels.native_order];
low_impedance_idx = all_impedances < impedance_threshold;
low_impedance_ch = native_order(low_impedance_idx);

Fs= Intan_struct.frequency_parameters.amplifier_sample_rate;
samplerate = Fs;
desired_length = ceil(length(concat_data)/Fs)*Fs;
nZEROS = desired_length - length(concat_data);
concat_data = concat_data(low_impedance_ch+1,:); % add 1 for indexing
concat_data_appended = [concat_data zeros(size(concat_data,1),nZEROS)];

writemda(concat_data_appended,fullfile(input_dir,'data.mda'),'int16');

clear concat data concat_data_appended amplifier_data

% .mda file name for Mountainsorrt 
files = dir(fullfile(input_dir,'*.mda')); file_names ={files.name};
file_name = file_names{1};

% Change geom.csv based on recorded channels
cd(working_dir)
[~,~,good_ch] = intersect(low_impedance_ch,sort(channel_map));
map = generateMap(channel_map,good_ch,input_dir);