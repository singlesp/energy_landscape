function [tri cap_edge] = HeadMeshTriangles(elec_pos, min_angle, max_dist_rel)
% [tri cap_edge] = HeadMeshTriangles(elec_pos, min_angle, max_dist_rel)
%
% elec_pos = Nx3 electrode positions
% min_angle = remove any questionable triangles have angles below this
%   (default 30 degrees)
% max_dist_rel = remove any questionable triangles with at least two sides
%   that are >= max_dist_rel*(average electrode neighbor spacing)
%   (default 2)
%
% tri = triangle verts connecting electrodes
% cap_edge = indices of electrodes along the edge of the cap

if(nargin < 2 || isempty(min_angle))
    min_angle = 30;
end

if(nargin < 3 || isempty(max_dist_rel))
    max_dist_rel = 2;
end

num_elec = size(elec_pos,1);
% 
% x = chanpos(:,1);
% y = chanpos(:,2);
% z = chanpos(:,3);
% 
% [th,ph,r] = cart2sph(x,y,z); % "project" points onto sphere
% 
% %repeat points at theta+2pi to force mesh to wrap
% th2 = [th; th+2*pi]; 
% ph2 = [ph; ph];
% 
% tri = delaunay(th2,ph2); %form triangle mesh
% 
% %adjust triangle vertex indices to refer to ORIGINAL points (not copied/shifted)
% tri(tri > numchan) = tri(tri > numchan) - numchan;
% 
% %sort each triangle's vertex indices and throw out duplicates
% tri = sort(tri,2); 
% tri = unique(tri,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = distance(elec_pos',elec_pos');
Ds = sort(D,2);
d = mean(Ds(:,2)); %average distance to nearest neighbor
avg_nei_dist = d;

%%%%%% collapse onto sphere
[c r] = spherefit(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3));
elec_sph = elec_pos - repmat(c',num_elec,1);
elec_sph = normrows(elec_sph);

D = distance(elec_sph',elec_sph');
Ds = sort(D,2);
d = mean(Ds(:,2)); %average distance to nearest neighbor (on sphere)
avg_nei_dist_sph = d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% load spherical surface to find sphere points near cap edge
sph_file = 'sphere_surf_simp.mat';
Sph = load(sph_file);
sph_loc = Sph.POS;
sph_tri = Sph.TRI1;

[c r] = spherefit(sph_loc(:,1),sph_loc(:,2),sph_loc(:,3));
sph_loc = sph_loc - repmat(c',size(sph_loc,1),1);
sph_loc = normrows(sph_loc);

sph_loc = sph_loc * .95; %shrink it a little
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

khull = convhulln(elec_sph);

%find scalp triangles that fall within this hull
tri_incap = inhull(sph_loc(sph_tri(:),:),elec_sph,khull);
tri_incap = reshape(tri_incap,size(sph_tri));
tri_incap = all(tri_incap,2);

%find the open edges of this set of faces (the edge of the cap scalp surface)
openedge = surfedge(sph_tri(tri_incap,:));
loops = extractloops(openedge); %breaks into connected loops

loop_breaks = DivideNanVector(loops);
num_loops = size(loop_breaks,1);

sph_edge = [];
for i = 1:num_loops
    loop = loops(loop_breaks(i,1):loop_breaks(i,2));
    sph_edge = [sph_edge; sph_loc(loop,:)];
end

%find the closest cap verts to each surface (sphere) edge vert
D = distance(sph_edge',elec_sph');
[Ds idx] = sort(D,2);
cap_edge = unique(idx(:,1));

%%%%%% find triangle verts that are in cap_edge
idx = ismember(khull(:),cap_edge);

%%%%%% find triangles made up of 3 cap edge verts
outside_tri = false(size(khull));
outside_tri(idx) = true;
outside_tri = all(outside_tri,2);
edge_tri = khull;

%%%%%% calculate all triangle edge distances and angles
a = rownorm(elec_sph(edge_tri(:,1),:)-elec_sph(edge_tri(:,2),:));
b = rownorm(elec_sph(edge_tri(:,1),:)-elec_sph(edge_tri(:,3),:));
c = rownorm(elec_sph(edge_tri(:,2),:)-elec_sph(edge_tri(:,3),:));

A = (a.^2-b.^2-c.^2) ./ (-2*b.*c);
B = (b.^2-a.^2-c.^2) ./ (-2*a.*c);
C = (c.^2-a.^2-b.^2) ./ (-2*a.*b);

%all triangles with an angle < 20deg
skinny_tri = any([A B C] > cos(min_angle*pi/180),2); 

%all triangles where at least 2 sides are > twice the average neighbor dist
dthresh = max_dist_rel*avg_nei_dist_sph;
large_tri = sum([a b c] > dthresh,2) > 1;

%all triangles made up of 3 cap_edge  verts that are either skinny or large
outside_tri = outside_tri & (skinny_tri | large_tri); %remove 

khull = khull(~outside_tri,:);

tri = khull;

%now find the new cap_edge defined by the edges of the new triangles
cap_edge = surfedge(khull);
cap_edge = unique(cap_edge(:));