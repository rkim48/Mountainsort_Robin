function [amplifier_data,amplifier_channels,frequency_parameters,ADC] = load_IntanRHS_data_wrapper(filename,selection)

 
[amplifier_data,amplifier_channels,frequency_parameters,ADC]=read_Intan_RHS2000_file(filename,selection);

end