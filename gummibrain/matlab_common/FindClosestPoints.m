function [orig_idx match_err] = FindClosestPoints(new_verts, orig_verts, varargin)
% [orig_idx match_err] = FindClosestPoints(new_verts, orig_verts, [status_update = false])
%
% for each point in new_verts, return the index of the closest point in
% orig_verts

status_update = false;
if(nargin > 2)
    status_update = varargin{1};
end

num_verts = size(new_verts,1);

orig_idx = zeros(num_verts,1);
match_err = zeros(num_verts,1);

match_err_thresh = 0.05; %not currently used

for i = 1:num_verts
    new_to_old_dist = sum([orig_verts(:,1)-new_verts(i,1) orig_verts(:,2)-new_verts(i,2) orig_verts(:,3)-new_verts(i,3)].^2,2);
    [d idx] = min(new_to_old_dist);

    orig_idx(i) = idx;
    match_err(i) = sqrt(d);
    
    if(status_update && ~mod(i,500))
        fprintf('%.0f%%... ',100*i/num_verts);
    end
    
    if(status_update && ~mod(i,2500))
        fprintf('\n');
    end
end

if(status_update)
    fprintf('Done!\n');
    fprintf('Average match error: %g\n', mean(match_err));
    fprintf('Largest match error: %g\n', max(match_err));
end