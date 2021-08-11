function result = MaskSphereSurface(elecpos,elecmask,surf)
% give electrode locations and surface
% return surface mask and mask_edge
% assumes elecpos already aligned to surf

[c, r] = spherefit(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% 3d volume from center of sphere to electrodes in the mask 
%   (extended out to make sure they are outside sphere)
pc = [elecpos(elecmask,:)*2; c.'];
fc = convhull(pc);
mask2 = inhull(surf.vertices, pc,fc).';

%3d distance between verts
D = distance(elecpos',elecpos');
Ds = sort(D,2);
avg_nei_dist = median(Ds(:,2));

%convert spatial dist to surface distance (on unit sphere)
D = 2*r*asin(D/(2*r)); 
Ds = sort(D,2);
avg_nei_dist_sph = median(Ds(:,2));

%extend mask by half of electrode spacing
D = distance(surf.vertices(mask2,:).',surf.vertices.');
%sph_capmask = min(Da,[],1);
%mask2 = +(min(D,[],1) < .5*avg_nei_dist_sph);
mask2 = +(min(D,[],1) < .5*avg_nei_dist);

nei = vertex_neighbours(surf);
mask2_exp = mesh_diffuse(mask2,nei,1)>0;
cap_faces = sum(mask2_exp(surf.faces),2)==3;

openedge = surfedge(surf.faces(cap_faces,:));
loops = extractloops(openedge); %breaks into connected loops

loop_breaks = DivideNanVector(loops);
num_loops = size(loop_breaks,1);

mask_edge = [];
for i = 1:num_loops
    loop = loops(loop_breaks(i,1):loop_breaks(i,2));
    mask_edge = [mask_edge; surf.vertices(loop([1:end 1]),:)];
end

result = struct('mask',mask2','mask_edge',mask_edge);