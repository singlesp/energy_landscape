function trimfile(filename,N)

M = [];
fid = fopen(filename,'r');
tline = fgetl(fid);
while(ischar(tline))
    M{end+1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

if(N < 0)
    M = M(1:end+N);
elseif(N < numel(M))
    M = M(1:N);
end

fid = fopen(filename,'w');
for n = 1:numel(M)
    fprintf(fid,'%s\n',M{n});
end
fclose(fid);
