function M = MountainsortParams
% MountainsortParams Specify Mountainsort parameters and pass values to external bash script that runs Mountainsort
%   Returns a structure M containing values

    working_dir =  '/media/robin/Shared/mountainsort_Robin'; 
    input_dir = uigetdir('server','Select folder with data');
    % If you read data from server, Mountainsort can not write outputs to
    % server so instead write to a local directory 
    output_dir = uigetdir('','Select folder to save outputs to');
    file_name = 'data.mda';
    samplerate=20000;
    geom_file = 'geom.csv';
    clip_size=60;
    freq_min=300;
    freq_max=5000;
    detect_interval=27;
    detect_threshold=4.5;
    adjacency_radius=100;
    detect_sign=-1;
    
    % Pass values to copy of Mountainsort bash script 
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
    
    % return structure with values to use in other functions 
    M.working_dir =  working_dir;
    M.input_dir = input_dir;
    M.output_dir = output_dir;
    M.file_name = file_name;
    M.samplerate = samplerate;
    M.geom_file = geom_file;
    M.clip_size = clip_size;
    M.freq_min = freq_min;
    M.freq_max = freq_max;
    M.detect_interval = detect_interval;
    M.detect_threshold = detect_threshold;
    M.adjacency_radius = adjacency_radius;
    M.detect_sign = detect_sign;
end