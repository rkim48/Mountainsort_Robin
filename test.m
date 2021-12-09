addpath('mdaio/')
%% Use conda m-file to switch to conda environment with Mountainsort code

% Initialize
conda.init
% List available conda environments
conda.getenv
% Set environment
conda.setenv('mountainlab2')

%% Change ML config variables if needed 

% system('ml-config')
% system('ML_PACKAGE_SEARCH_DIRECTORY=/home/robin/miniconda3/envs/mountainlab2/etc/mountainlab/packages')
% system('ML_CONFIG_FILE=/home/robin/miniconda3/envs/mountainlab2/etc/mountainlab/mountainlab.env')

%%
working_dir = pwd;
% input_dir = uigetdir;
output_dir = fullfile(input_dir,'outputs'); mkdir(output_dir);
files = dir(fullfile(input_dir,'*.mda')); file_names ={files.name};
file_name = file_names{1};
samplerate=20000;
geom_file = 'geom.csv';
clip_size=60;
freq_min=300;
freq_max=5000;
detect_interval=27;
detect_threshold=4.5;
adjacency_radius=100;
detect_sign=-1;

bash_script = 'Mountainsort_Robin.sh';
execute_dir = '~/miniconda3/envs/mountainlab2/work_folder';
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

% Write changes to .sh file
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

cd(execute_dir);
system(strcat('./',copy_str));

cd(working_dir);