function img = mri_orthmip(vol,dims,bb,do_minimum)
if(nargin < 2 || isempty(dims))
    dims = 1:3;
end
if(nargin < 3)
    bb = [];
end
if(nargin < 4)
    do_minimum = false;
end
if(ischar(do_minimum))
    if(strncmp('minimum',do_minimum,numel(do_minimum)))
        do_minimum = true;
    else
        do_minimum = false;
    end
end


mipdim = [3 2 1];
%mipdim = [1 2 3];
M = mri_orthvol(vol,bb);
if(do_minimum)
    M = -M;
end

img = {};
for d = 1:numel(dims)
    posmax = squeeze(max(M,[],mipdim(dims(d))));
    negmax = squeeze(min(M,[],mipdim(dims(d))));
    posmax(abs(negmax) > abs(posmax)) = negmax(abs(negmax) > abs(posmax));
    img{d} = posmax;
    if(dims(d) == 3)
        img{d} = flipud(img{d});
    end
    if(do_minimum)
        img{d} = -img{d};
    end
end


% return;
% %%
% if(nargin < 3 || isempty(bb))
%     bb = .5*[-256 -256 -256; 256 256 256];
% end
% 
% if(ischar(vol))
%     vol = spm_vol(vol);
% end
% 
% 
% img = {};
% for d = 1:numel(dims)
%     dim = dims(d);
%     spmvec = [0 0 0];
%     spmvec(dim) = 1;
%     mipimg = [];
%     for i = floor(bb(1,dim)):ceil(bb(2,dim))
%         tmp = mri_orthslice(vol,spmvec*i,dim,bb);
%         if(isempty(mipimg))
%             mipimg = tmp{1};
%         else
%             mipimg = max(mipimg,tmp{1});
%         end
%     end
%     img{d} = mipimg;
% end
