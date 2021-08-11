function [cx cy cz] = compute_mesh_contour(tri,verts,vals,level)
% [cx cy cz] = mesh_contour(tri,verts,vals,level)
% 
% returns a set of connected line segments representing the contour of a triangle mesh
% at the given contour level
%
% tri: Mx3 set of vert indices
% verts: Nx3 set of [x y z]
% vals: Nx1 set of values for each vertex
% level: level at which to find contour
%
% cx,cy,cz = Px1 containing x,y,z a set of disconnected paths making up the contour.  
%   if a contour has multiple disconnected paths for a given level, these
%   paths are listed sequentially, divided by a NaN
% 
% Segments can plot as follows:
%     [cx cy cz] = mesh_contour(tri,verts,vals,level)
%     path_breaks = DivideNanVector(cx);
%     for i = 1:size(path_breaks,1)
%         p = path_breaks(i,1):path_breaks(i,2);
%         plot3(cx(p),cy(p),cz(p));
%     end

%create edge list such that a triangle's edges are grouped
%and fix(edgeindex/3)+1 gives you the triangle
edges = zeros(3*size(tri,1),2);
edges(1:3:end,:) = tri(:,[1 2]);
edges(2:3:end,:) = tri(:,[2 3]);
edges(3:3:end,:) = tri(:,[1 3]);

%%% ed_u = ed(m,:),   ed = ed_u(n)
%[ed_u m n] = unique(sort(edges,2),'rows');

%cross_edge_idx = (vals(edges(:,1)) > thresh) & (vals(edges(:,2)) < thresh);
cross_edge_idx = find(sum(vals(edges) > level,2) == 1);
%cross_edge_tri_idx = edge_tri_idx(cross_edge_idx,:);

%also find edges where one vertex is an exact match
%and figure out a way to deal with those connections
%cross_edge_idx = sum(c(edges) == thresh,2) == 1;

if(isempty(cross_edge_idx))
    cx = [];
    cy = [];
    cz = [];
    return;
end

cross_vert_idx = edges(cross_edge_idx,:);

cross_vert_val = vals(cross_vert_idx);
dist_to_thresh = abs(cross_vert_val - level);
w_thresh = 1 - dist_to_thresh ./ repmat(sum(dist_to_thresh,2),1,2);
w_thresh(dist_to_thresh < 1e-13) = 1;

cross_x = [verts(cross_vert_idx(:,1),1) verts(cross_vert_idx(:,2),1)];
cross_y = [verts(cross_vert_idx(:,1),2) verts(cross_vert_idx(:,2),2)];
cross_z = [verts(cross_vert_idx(:,1),3) verts(cross_vert_idx(:,2),3)];

cross_x = sum(cross_x.*w_thresh,2);
cross_y = sum(cross_y.*w_thresh,2);
cross_z = sum(cross_z.*w_thresh,2);

idx = reshape(1:numel(cross_edge_idx),2,numel(cross_edge_idx)/2)';
% cx = cross_x(idx);
% cy = cross_y(idx);
% cz = cross_z(idx);
% return;

%%% if an edge shows up multiple times (eg, in adjacent triangles), 
%%% keep it but replace the indices so they all refer to the same
%%% edge.  (ex: if edge 1 and 4 are the same, replace all 4's with 1's)
[cvidx_u m n] = unique(sort(cross_vert_idx,2),'rows'); 
idx_u = m(n(idx));

paths = extractpaths(idx_u);
paths_nans = isnan(paths);

paths_denan = paths;
paths_denan(paths_nans) = 1;

cx = cross_x(paths_denan);
cy = cross_y(paths_denan);
cz = cross_z(paths_denan);

cx(paths_nans) = nan;
cy(paths_nans) = nan;
cz(paths_nans) = nan;

