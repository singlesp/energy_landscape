function img = mri_orthvol(vol,bb)


if(ischar(vol))
    vol = spm_vol(vol);
end

if(nargin < 2 || isempty(bb))
    bb = spm_get_bbox(vol,'fv');
    %bb = .5*[-256 -256 -256; 256 256 256];
end

img = [];

dim = 1;
spmvec = [0 0 0];
spmvec(dim) = 1;

i1 = floor(bb(1,dim));
i2 = ceil(bb(2,dim));
for i = i1:i2
    tmp = mri_orthslice(vol,spmvec*i,dim,bb);
    if(isempty(img))
        img = zeros([size(tmp{1}) i2-i1+1]);
    end
    img(:,:,i-i1+1) = tmp{1};
end
