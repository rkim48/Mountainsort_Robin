% Intan channel 38 is the shallowest electrode while 25 is the deepest. 
% Let channel 38 be at position (0,0). 
function [map, depth_idx] = generateMap(channel_map,good_ch,save_dir)
%     rhd_channels  = [25 47 27 45 29 43 31 41 23 33 21 35 19 37 17 39 24 46 26 44 28 42 30 40 22 32 20 34 18 36 16 38];
%     rhs_channels = [27 16 26 17 25 18 24 19 28 23 29 22 30 21 31 20 4 15 5 14 6 13  7  12  3  8  2  9  1 10  0 11];
    nCH = numel(channel_map);
    pitch = 60; % um
    x = zeros(nCH,1);
    y_presort =zeros(nCH,1);
    for i = 1:nCH
        y_presort(i) =pitch * (i-1); 
    end
    
    pos_arr = [(1:nCH)' flip(channel_map) y_presort]; % indices, Intan channels, y_coordinates

    [sorted_Intan_ch,sorted_idx] = sort(channel_map);
    y = y_presort(sorted_idx);
    depth_idx = (y / pitch) + 1;
    
    ordered_pos_arr =  [(1:nCH)' sorted_Intan_ch y];  % pos_arr but ordered with respect to Intan channel order

    map = [x y];
    map = map(good_ch,:);
    depth_idx = depth_idx(good_ch);
    
    writematrix(map,fullfile(save_dir,'geom.csv'))
end