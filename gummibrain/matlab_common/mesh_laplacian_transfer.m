function M_laplace = mesh_laplacian_transfer(elec_pos, elec_tri, h, h_type)
% M_laplace = mesh_laplacian_transfer(elec_pos, elec_tri, h, h_type)
%
% Note: can also be used on splined surfaces
%
% Using the mesh laplacian operator from:
% Belkin, M., Sun, J. & Wang, Y. "Discrete laplace operator on meshed surfaces."
% Proceedings of the twenty-fourth annual symposium on Computational geometry 278-287 
% (2008).doi:10.1145/1377676.1377725

num_elec = size(elec_pos,1);

[c r] = spherefit(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3));
elec_sph = elec_pos - repmat(c',num_elec,1);
elec_sph = normrows(elec_sph);

a = rownorm(elec_sph(elec_tri(:,1),:)-elec_sph(elec_tri(:,2),:));
b = rownorm(elec_sph(elec_tri(:,1),:)-elec_sph(elec_tri(:,3),:));
c = rownorm(elec_sph(elec_tri(:,2),:)-elec_sph(elec_tri(:,3),:));

sp = sum([a b c],2)/2;
tri_areas = sqrt(sp.*(sp-a).*(sp-b).*(sp-c));



D = distance(elec_sph',elec_sph'); %3d distance between verts
Ds = sort(D,2);
avg_nei_dist = mean(Ds(:,2));
    
r=1;
D = 2*r*asin(D/(2*r)); %convert spatial dist to surface distance (on unit sphere)
Ds = sort(D,2);
avg_nei_dist_sph = mean(Ds(:,2));
    

if(nargin < 4 || isempty(h_type))
    h_type = 'relative';
end

if(nargin < 3 || isempty(h))
    h = avg_nei_dist_sph/10;
    h_type = 'absolute';
end

if(strcmpi(h_type,'relative'))
    h = h*avg_nei_dist_sph;
end

areamat = [tri_areas; tri_areas; tri_areas];

distmat = exp(-D.^2/(4*h)); %convert to weight used in formula

trivert = elec_tri(:);
trivert_dist = distmat(elec_tri,:);

M_laplace = zeros(num_elec,num_elec);

for w=1:num_elec
    A = -trivert_dist;
    A(trivert ~= w,:) = 0;
    A(:,w) = (trivert ~= w).*trivert_dist(:,w);
    
    M_laplace(:,w) = A'*areamat;
end
M_laplace = M_laplace/(4*3*pi*h^2);