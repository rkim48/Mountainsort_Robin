
ChMapNum = 2
load('Dat_V_Map.mat')

% load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
switch ChMapNum
   case 1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 0
   load('ch_map_pink.mat');
   case 2
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_large.mat');
   case 3
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_Small.mat');
    case 4
    load('D:\Box Sync\0919 stroke\Finger_map.mat');
    case 5
    load('Mirro_Oversampling_hippo_map.mat');
end

Ch_Map=Ch_Map_new;
if ChMapNum==5 || ChMapNum==6 || ChMapNum==7 || ChMapNum==8 
Ch_Map=Maps{sfol};
end

if ChMapNum==9
Ch_Map=[Maps{sfol}(:,1) ;Maps{sfol}(:,2) ];
end


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
% Ox=12.5;    
% Oy=12.5;
%         gap=5;
%         Gw=25;
%         Gh=25;
 if Hip~=1       
    Ox=15;    
    Oy=15;
    gap=20;
    Gw=30;
    Gh=30;

    numGx = 8;
    numGy = 4;
    Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
    Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
    [X,Y] = meshgrid(Gx,Gy);
        
    CSVList = zeros(size(Dat_V_Map));
    for ele = 1:size(Dat_V_Map,1)
        [row,col]=find(Ch_Map_2==ele);
        if ~isempty(row)
            CSVList(ele,1)=X(row,col);
            CSVList(ele,2)=Y(row,col);
        end
    end
    if mode ==1
        csvwrite([global_save_path str '/' 'geom.csv'],CSVList);
    else
        csvwrite(['geom.csv'],CSVList);    
    end
 end
      
      %%
      if Hip==1
Ox=10;    
Oy=10;
gap=10;
Gw=20;
Gh=20;
        if ChMapNum==7 || ChMapNum==8
Ox=10;    
Oy=10;
gap=5;
Gw=15;
Gh=15;
end

numGx = 2;
numGy = 16;

if ChMapNum==9
Ox=10;    
Oy=10;
gap=45;
Gw=15;
Gh=15;
numGx = 32;
numGy = 1;
end




Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);
 if size(X,1)==1
     X=X';
     Y=Y';
 end
      CSVList = zeros(size(Dat_V_Map));
      for ele = 1:size(Dat_V_Map,1)
      [row,col]=find(Ch_Map_2==ele);
      CSVList(ele,1)=X(row,col);
      CSVList(ele,2)=Y(row,col);
      end
      
      if mode ==1
      csvwrite([global_save_path str '/' 'geom.csv'],CSVList);
      else
      csvwrite(['geom.csv'],CSVList);    
      end      
      
      end
      
      