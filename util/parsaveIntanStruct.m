function parsaveIntanStruct(amplifier_channels,frequency_parameters,working_dir)
    
    Intan_struct.amplifier_channels = amplifier_channels;
    Intan_struct.frequency_parameters = frequency_parameters;
    save(fullfile(working_dir,'Intan_struct'),'Intan_struct');

end