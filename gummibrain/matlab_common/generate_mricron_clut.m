%cmap = GenerateColormap([0 .4 .5 .6 1],[0 1 0; 0 0 1; 0 0 0; 1 1 0; 1 0 0],256);
cmap = GenerateColormap([0 .5 1],[0 1 0; 0 0 0; 1 0 0],256);
X = [[0:size(cmap,1)-1]' round(cmap*255)];

clut_name = 'plusminus4';

clut_file = sprintf('~/Applications/mricron/lut/%s.lut',clut_name);
fid = fopen(clut_file,'w');
fprintf(fid,'* s=byte\tindex\tred\tgreen\tblue\n');
fprintf(fid,'S\t%d\t%d\t%d\t%d\n',X');
fclose(fid);