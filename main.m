% Main script 
addpath('mdaio/')
addpath('channel_maps/')
addpath('util/')
%% Run Mountainsort
cd(working_dir);
system(strcat('./',copy_str));
%%
cd(working_dir);
animalID = 'ICMS38'; % for figure title
output_metrics_name = extractClusterMetrics(output_dir,ch_idx_table,file_name,map,Fs);
output_curated_name = curateClusters(output_metrics_name,output_dir,Fs);
visualizeClusters(output_curated_name,output_dir,Fs,animalID);


% manualCuration(output_curated_name,output_dir,Fs,animalID);
%%
cd(execute_dir)
cd(output_dir)