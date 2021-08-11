function M = idx2vol(volsize,idx,vals,emptyval)
% M = idx2vol(volsize,idx,vals,emptyval)

if(nargin < 3 || isempty(vals))
    vals = 1;
end

if(numel(vals) == 1)
    vals = vals*ones(numel(idx),1);
end

if(nargin < 4 || isempty(emptyval))
    emptyval = 0;
end

if(max(idx) > prod(volsize))
    error('idx cannot exceed volume size');
end

if(numel(idx) == numel(vals))
    M = emptyval*ones(volsize);
    M(idx) = vals;
elseif(numel(idx) == size(vals,1))
    valdim=size(vals);
    nvox=prod(volsize);
    M = emptyval*ones([nvox valdim(2:end)]);
    M(idx,:)=reshape(vals,numel(idx),[]);
    M=reshape(M,[volsize valdim(2:end)]);
else
    error('vals and idx must be the same size');
end


