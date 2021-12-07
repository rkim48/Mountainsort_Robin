fid = fopen('sftest.py','wt');
str0='from spikesorters import mountainsort4';
fprintf(fid,str0);
fprintf(fid,'\n');

str1='mountainsort4.execute(';
fprintf(fid,str1);
fprintf(fid,'\n');

str2=["recording_dir='" pwd "',"];
str2=strcat(str2{:});
fprintf(fid,str2);
fprintf(fid,'\n');

str3=["firings_out='" pwd "/firings.mda',"];
str3=strcat(str3{:});
fprintf(fid,str3);
fprintf(fid,'\n');

str4="detect_sign=0,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="whiten=False,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="clip_size=54,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="detect_interval=27,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="freq_min=0,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="freq_max=0,";
fprintf(fid,str4);
fprintf(fid,'\n');


str4="detect_threshold=45,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="adjacency_radius=100,";
fprintf(fid,str4);
fprintf(fid,'\n');

str4="_container=None)";
fprintf(fid,str4);
fprintf(fid,'\n');

fclose(fid);