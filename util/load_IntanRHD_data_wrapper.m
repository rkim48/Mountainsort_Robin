function [amplifier_data,amplifier_channels,frequency_parameters,channel_struct] = load_IntanRHD_data_wrapper(filename,selection)

 
[amplifier_data,amplifier_channels,frequency_parameters,channel_struct]=read_Intan_RHD2000_fileV3(filename,selection);

end