% Main script 
%% Run Mountainsort
cd(execute_dir);
system(strcat('./',copy_str));
%%
cd(working_dir);
output_metrics_name = extractClusterMetrics(output_dir,file_name,map,Fs);
output_curated_name = curateClusters(output_metrics_name,output_dir,Fs);
visualizeClusters(output_curated_name,output_dir,Fs);
%%
cd(execute_dir)
