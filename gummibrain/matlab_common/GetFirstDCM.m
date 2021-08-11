function dcmfile = GetFirstDCM(filename,scan)

fid = fopen(filename,'r');
M = textscan(fid,'%d %s');
fclose(fid);

dcmfile = M{2}{M{1} == scan};