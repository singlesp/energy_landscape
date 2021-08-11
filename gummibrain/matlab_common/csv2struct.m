function S = csv2struct(filename,fmt)

fid=fopen(filename);

hdrtxt=strsplit(fgetl(fid),',','collapsedelimiters',false);

if(~exist('fmt','var') || isempty(fmt))
    fmt=repmat('%s,',1,numel(hdrtxt));
    fmt=fmt(1:end-1);
end

M = {};
while(~feof(fid))
    line=strsplit(fgetl(fid),',','collapsedelimiters',false);
    if(isempty(line))
        continue
    end
    if(~isempty(M) && numel(line) > size(M,2))
        M = [M repmat({''},size(M,1),numel(line)-size(M,2))];
    end
    if(numel(line) < size(M,2))
        line = [line repmat({''},1,size(M,2)-numel(line))];
    end
    M = [M; line];
end
fclose(fid);

a=M;