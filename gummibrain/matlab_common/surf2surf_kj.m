function newvals = surf2surf_kj(surfvals, nn12, nn21, mode)

if(size(surfvals,1) == numel(nn12) && size(surfvals,2) == numel(nn12))
    surfvals = surfvals.';
end

if(nargin < 3)
    nn21 = [];
end

if(nargin < 4 || isempty(mode))
    mode = 'nnfr';
end

newvals = surfvals(nn12,:);

if(strcmpi(mode,'nnfr') && ~isempty(nn21))
    is_vert_used = false(size(surfvals,1),1);
    is_vert_used(nn12) = true;
    newvals(nn21(~is_vert_used),:) = (newvals(nn21(~is_vert_used),:) + surfvals(~is_vert_used,:))/2;
end