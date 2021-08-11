function outfile = convert_bvanalyzer2spmanalyzer(bvhdrfile, outfile)
if(nargin < 2)
    outfile = '';
end

[bvdir,bvfile,bvext] = fileparts(bvhdrfile);
if(strcmpi(bvext,'.img'))
    bvext = '.hdr';
end

if(isempty(outfile))
    outfile = sprintf('%s/%s_SPM%s',bvdir,bvfile,bvext);
else
    [newdir,newfile,newext] = fileparts(outfile);
    if(isempty(newext))
        newext = bvext;
    end
    if(strcmpi(newext,'.img'))
        newext = '.hdr';
    end
    if(isempty(newdir))
        newdir = bvdir;
    end

    if(isempty(newfile))
       newfile = sprintf('%s_SPM%s',bvfile,bvext);
    end
    outfile = sprintf('%s/%s%s',newdir,newfile,newext);
end

%%%%%%%%%%%%%%%%%%%%%%%%%

X = xff(bvhdrfile);

isvmr = false;
M = [];
try
    M = X.VoxelData;
catch
end

if(isempty(M))
    try
        M = X.VMRData;
        isvmr = true;
    catch
    end
end

pre_type = class(M);
M = permute(double(M),[3 1 2]);
M = M(:,end:-1:1,end:-1:1);
%M = flipdim(flipdim(flipdim(M,1),2),3);
M = cast(round(M),pre_type);
%X.ImgDim.PixSpacing(1) = -1;

if(isvmr)
    X.VMRData = M;
    X.ImgDim.Dim(2:4) = size(M);
else
    X.VoxelData = M;
    X.ImgDim.Dim(2:4) = size(M);
end

X.SaveAs(outfile);
X.ClearObject();
