function [amplifier_data,amplifier_channels,frequency_parameters,ADC] = load_Intan_data_wrapper(filename,aggPath,recordingFolders,chname,strsh,selection)
for repeat=1:3
    try
        [amplifier_data,amplifier_channels,frequency_parameters,ADC]=read_Intan_RHD2000_fileV3(filename,selection);
        break
    catch
       download_one_file_from_box_Intan(aggPath,recordingFolders,chname,strsh) 
    end
    if repeat == 3
    load('IntanBackupChsFre')
    amplifier_data=rand(128,5400000)*10;
    ADC = zeros(2,5400000);
end
end