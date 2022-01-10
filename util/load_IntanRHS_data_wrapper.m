function [amplifier_data,amplifier_channels,frequency_parameters,channel_struct] = load_IntanRHS_data_wrapper(filename)

 
[amplifier_data,amplifier_channels,frequency_parameters,channel_struct]=read_Intan_RHS2000_file(filename);

end