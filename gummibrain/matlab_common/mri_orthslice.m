function img = mri_orthslice(V,center,dims,bb)
if(nargin < 2 || isempty(center))
    center = [0 0 0];
end

if(nargin < 3 || isempty(dims))
    dims = 1:3;
end
if(nargin < 4 || isempty(bb))
    bb = .5*[-256 -256 -256; 256 256 256];
end

if(ischar(V))
    V = spm_vol(V);
end

sp = eye(4); %st.Space
premul = eye(4);
Dims = round(diff(bb)'+1);
interpmode = 1; %interp. 0=nearest neighbor

M = sp\premul*V.mat;

TM0 = [ 1 0 0 -bb(1,1)+1
    0 1 0 -bb(1,2)+1
    0 0 1 -center(3)
    0 0 0 1];
TM = inv(TM0*M);
TD = Dims([1 2]);

CM0 = [ 1 0 0 -bb(1,1)+1
    0 0 1 -bb(1,3)+1
    0 1 0 -center(2)
    0 0 0 1];
CM = inv(CM0*M);
CD = Dims([1 3]);

SM0 = [ 0 -1 0 +bb(2,2)+1
    0  0 1 -bb(1,3)+1
    1  0 0 -center(1)
    0  0 0 1];
SM = inv(SM0*M);
SD = Dims([2 3]);

img = {};
for d = 1:numel(dims)
    if(dims(d) == 1)
        img{d} = spm_slice_vol(V,SM,SD,interpmode)';
    elseif(dims(d) == 2)
        img{d} = spm_slice_vol(V,CM,CD,interpmode)';
    elseif(dims(d) == 3)
        img{d} = spm_slice_vol(V,TM,TD,interpmode)';
    end
end



