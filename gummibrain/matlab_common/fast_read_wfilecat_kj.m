function w = fast_read_wfilecat_kj(fname,numvertex)

if(nargin < 2 || isempty(numvertex))
    numvertex = 0;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open w file %s.', fname) ;
  error(str) ;
end

numcat = str2num(fgetl(fid));
if(~isempty(numcat))
    for i = 1:numcat+1
        tmp=fgetl(fid);
    end
end
wcat = {};
i = 0;
while(1)
    i = i+1;
    fread(fid, 1, 'int16') ;  % Skip ilat
    vnum = fast_fread3(fid) ; % Number of non-zero values
    if(isempty(vnum))
        break;
    end
    A = uint32(fread(fid,vnum*7,'uchar'));
    

    v = bitshift(A(1:7:end), 16) + bitshift(A(2:7:end),8) + A(3:7:end);
    w0 = bitshift(A(4:7:end),24) + bitshift(A(5:7:end),16) + bitshift(A(6:7:end),8) + A(7:7:end);

    numvertex = max(max(v),numvertex);
    wcat{i} = zeros(numvertex,1);
    wcat{i}(v+1) = typecast(w0,'single');
    
    if(feof(fid))
        break;
    end
end
fclose(fid);
wcat = cellfun(@(x)padto(x,numvertex),wcat,'uniformoutput',false);
w = cat(2,wcat{:});

end

function y = padto(x,len)
    sz = size(x);
    if(sz(1) == 1)
        y = zeros(1,len);
    else
        y = zeros(len,1);
    end
    y(1:numel(x)) = x;
end