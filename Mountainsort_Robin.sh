input_dir=/media/robin/Shared/ICMS/ICMS19/June29_35samples
output_dir=/media/robin/Shared/ICMS/ICMS19/June29_35samples
samplerate=30000
file_name=raw.mda
geom_file=map.csv
clip_size=54
freq_min=0
freq_max=0
detect_interval=27
detect_threshold=45
adjacency_radius=100
detect_sign=-1

# Run bandpass filter stage of Mountainsort
ml-run-process ephys.bandpass_filter \
	--inputs timeseries:$input_dir/$file_name \
	--outputs timeseries_out:$output_dir/filt_${file_name} \
	--parameters samplerate:$samplerate freq_min:$freq_min freq_max:$freq_max

ml-run-process ephys.whiten \
  --inputs timeseries:$output_dir/filt_${file_name} \
	--outputs timeseries_out:$output_dir/pre.mda.prv

# Spike sorting
# Specify the detect threshold in standard deviations
ml-run-process ms4alg.sort \
	--inputs \
		timeseries:$output_dir/pre.mda.prv geom:$input_dir/$geom_file \
	--outputs \
		firings_out:$output_dir/firings_${file_name} \
	--parameters \
		detect_sign:$detect_sign \
		adjacency_radius:$adjacency_radius \
        detect_interval:$detect_interval \
		detect_threshold:$detect_threshold \
		clip_size:$clip_size \
		num_features:15 \
		num_workers:24 \

# Compute cluster metrics
ml-run-process ephys.compute_cluster_metrics \
	--inputs \
		timeseries:$output_dir/pre.mda.prv firings:$output_dir/firings_${file_name} \
	--outputs \
		metrics_out:$output_dir/cluster_metrics.json \
	--parameters \
		samplerate:$samplerate \
		clip_size:$clip_size \
		refrac_msec:1 \

# Compute templates from just filtered data (preferred since we get amplitudes in uV)
ml-run-process ephys.compute_templates \
	--inputs \
		timeseries:$output_dir/filt_${file_name} firings:$output_dir/firings_${file_name} \
	--outputs \
		templates_out:$output_dir/templates_${file_name} \
	--parameters \
		clip_size:$clip_size

# Some cluster metrics
ml-run-process ms3.isolation_metrics \
	--inputs \
		timeseries:$output_dir/pre.mda.prv firings:$output_dir/firings_${file_name} \
	--outputs \
		metrics_out:$output_dir/isolation_metrics_out.json \
		pair_metrics_out:$output_dir/pair_metrics_out.json \
	--parameters \
		compute_bursting_parents:true

#combine metrics
ml-run-process ms3.combine_cluster_metrics \
	--inputs \
	metrics_list:$output_dir/cluster_metrics.json \
	metrics_list:$output_dir/isolation_metrics_out.json \
	--outputs \
	metrics_out:$output_dir/combine_metrics_new.json \

# Mountain View
# You can do manual curation of clusters using Mountain View
# Make sure to output the curated firings file in the GUI

# Change input data to filt.mda if you want to see clusters in uV instead of standard deviations
qt-mountainview --pre=$output_dir/pre.mda.prv \
		--firings=$output_dir/firings_${file_name} \
		--samplerate=$samplerate \
		--cluster_metrics=$output_dir/combine_metrics_new.json
