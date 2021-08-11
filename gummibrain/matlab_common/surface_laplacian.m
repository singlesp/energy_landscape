function [lap_vals lap_center_idx lap_surround_idx lap_shape lap_radius] = surface_laplacian(elec_vals, elec_pos, surf_vals, surf_verts, surf_norms, lap_shape, lap_radius, eval_pos)

if(nargin < 8)
    eval_pos = [];
end

if(nargin < 7 || isempty(lap_radius) || lap_radius < 0)
    D = distance(elec_pos',elec_pos');
    Ds = sort(D,2);
    d = mean(Ds(:,2));
    lap_radius = d; %for return val
else
    d = lap_radius;
end

elec_pos = [elec_pos; eval_pos];

num_elec = size(elec_pos,1);

%project electrode positions onto unit sphere
[c r] = spherefit(elec_pos);
% elec_pos = elec_pos - repmat(c',num_elec,1);
% elec_pos = normrows(elec_pos);

%find the closest surface vert for each electrode
D = distance(elec_pos',surf_verts');
[~, idx] = min(D,[],2);
elec_verts = surf_verts(idx,:);
elec_norms = surf_norms(idx,:);

lap_center_idx = idx;

if(mean(sum(elec_norms.*elec_pos,2)) < 0)
    %correct inward-pointing norms
    elec_norms = -elec_norms;
end

if(nargin < 6 || isempty(lap_shape))
    C = 20;
    lap_shape = C;  %for return val
else
    C = lap_shape; %how many points in circle
end



%r = 1;
alpha = 2*asin(d/(2*r));
cone_rad = r*sin(alpha);
%cone_height = r*cos(alpha);
cone_height = 1;  

th = linspace(0,2*pi,C+1);
th = th(1:C);
circ_pos = [zeros(1,C); cone_rad*cos(th); cone_rad*sin(th); ones(1,C)];
    
%generate 3 orthogonal axes at each electrode, where the first is surface
%normal
n1 = elec_norms;
n1 = normrows(n1);

%generate new non-parallel/degenerate point
n_tmp = n1; 
[~, i] = min(n_tmp,[],2); 
ind = sub2ind(size(n_tmp),1:size(n_tmp,1),i')';
n_tmp(ind) = n_tmp(ind)+max(n_tmp,[],2)+1; 
n2 = cross(n1,n_tmp,2);
n2 = normrows(n2);

%generate third norm by crossing first two
n3 = cross(n1,n2,2);
n3 = normrows(n3);

norm_mat = zeros(3*num_elec,4);
norm_mat(:,1) = reshape(n1',3*num_elec,1);
norm_mat(:,2) = reshape(n2',3*num_elec,1);
norm_mat(:,3) = reshape(n3',3*num_elec,1);
norm_mat(:,4) = reshape(elec_verts'*cone_height,3*num_elec,1); %shift after rotating

nc = norm_mat*circ_pos; %create laplacian surround copies

%reshape and reshuffle dimensions to accomodate next step
nc = reshape(nc,3,num_elec,C);
nc = permute(nc,[3 2 1]);
nc = reshape(nc,C*num_elec,3); %each 
% 
% 
%  figure;
% plot3(surf_verts(:,1),surf_verts(:,2),surf_verts(:,3),'.','markersize',1);
% hold on;
% plot3(elec_verts(:,1),elec_verts(:,2),elec_verts(:,3),'r.','markersize',10);
% axis vis3d equal;
% % 
% % size(nc)
%  nc2=nc+.01*rand(size(nc));
%  nc2 = nc;
%  plot3(nc2(:,1),nc2(:,2),nc2(:,3),'k-','markersize',1);

%find closest surface vert for each laplacian surround point for each
%electrode
D = distance(nc',surf_verts');
[~, lap_idx] = min(D,[],2);
lap_idx = reshape(lap_idx,C,num_elec);

if(isempty(elec_vals) || isempty(surf_vals))
    lap_vals = [];
else
    %lap_vert = vert(lap_vert_idx(:),:);
    lap_vert_vals = surf_vals(lap_idx);
    lap_vals = elec_vals - mean(lap_vert_vals,1);
end

if(nargout > 1)
    lap_surround_idx = lap_idx;
end