load('Dat_V_Map.mat')

switch ChMapNum
   case 1
       load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 2
       load('pavlo_1x32_intanMap_v2.mat');
end

Ch_Map=Ch_Map_new;

for i=1:size(Ch_Map,1)
    for j=1:size(Ch_Map,2)
        if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
    Ch_Map(i,j)=0;
    end
    end
end

Ch_Map_2 = Ch_Map; 
for i=1:size(Ch_Map,1)
    for j=1:size(Ch_Map,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
    Ch_Map_2(i,j)=0;
    else 
    Ch_Map_2(i,j)=find(Dat_V_Map(:,2)==Ch_Map(i,j));
    end
    end
end

% larger number is shallower 
spacing = 60;
y_value = [0:31] * spacing;
CSVList = zeros(size(Dat_V_Map));
for elec = 1:size(Dat_V_Map,1)
    [row,col]=find(Ch_Map_2==elec);
    if ~isempty(row)
        CSVList(elec,2)=y_value(elec);
    end
end

csvwrite(['geom.csv'],CSVList);    
 