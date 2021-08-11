function count = filelinecount(filename)

fid = fopen(filename,'rt');
count = 0;
while (ischar(fgetl(fid)))
  count = count + 1;
end
fclose(fid);
