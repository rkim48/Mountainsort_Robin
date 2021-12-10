% Intan channel 38 is the shallowest electrode while 25 is the deepest. 
% Let channel 38 be at position (0,0). 
function map = generateMap(Intan_ch,good_ch,save_dir)
%     Intan_ch = [25 47 27 45 29 43 31 41 23 33 21 35 19 37 17 39 24 46 26 44 28 42 30 40 22 32 20 34 18 36 16 38];
    Intan_ch = Intan_ch';
    nCH = numel(Intan_ch);
    pitch = 60; % um
    x = zeros(nCH,1);
    y_presort =zeros(nCH,1);
    for i = 1:nCH
        y_presort(i) =pitch * (i-1); 
    end
    pos_arr = [(1:nCH)' flip(Intan_ch') y_presort]; % indices, Intan channels, y_coordinates

    [sorted_Intan_ch,sorted_idx] = sort(flip(Intan_ch));
    y = y_presort(sorted_idx);

    ordered_pos_arr =  [(1:nCH)' sorted_Intan_ch' y];  % pos_arr but ordered with respect to Intan channel order

    map = [x y];
    map = map(good_ch,:);

    writematrix(map,fullfile(save_dir,'geom.csv'))
end