function w = fast_read_wfile_kj(fname,numvertex)

if(nargin < 2 || isempty(numvertex))
    numvertex = 0;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open w file %s.', fname) ;
  error(str) ;
end

fread(fid, 1, 'int16') ;  % Skip ilat
vnum = fast_fread3(fid) ; % Number of non-zero values

A = fread(fid,vnum*7,'uchar');
fclose(fid);

% % A needs to be a uint32 for these bitshifts to work properly
% A = uint32(A);
% v = bitshift(A(1:7:end), 16) + bitshift(A(2:7:end),8) + A(3:7:end);
% w0 = bitshift(A(4:7:end),24) + bitshift(A(5:7:end),16) + bitshift(A(6:7:end),8) + A(7:7:end);
% 
% numvertex = max(max(v),numvertex);
% w = zeros(numvertex,1);
% w(v+1) = typecast(w0,'single');

% A needs to be a uint8
A = uint8(A);
A = reshape(A,7,[]);
v = [A([3 2 1],:); zeros(1,size(A,2))];
v = typecast(v(:),'uint32');

w0 = A([7:-1:4],:);
w0 = typecast(w0(:),'single');
w = zeros(numvertex,1);
w(v+1) = w0;