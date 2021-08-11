function plotdipole(pos,vec,varargin)
% dipole([x y z],[nx ny nz])
% 
% %reorient if argument was [x y z]'
% if(size(pos,1) > size(pos,2))
%     pos = pos';
% end
% 
% if(size(vec,1) > size(vec,2))
%     vec = vec';
% end

size(pos)
size(vec)
scale_factor = 5;
if(nargin >= 7)
    scale_factor = varargin{7};
end

n = scale_factor*vec';

cyl_len = sqrt(n'*n);
sph_rad = cyl_len/3;
cyl_rad = cyl_len/6;
cyl_shift = sqrt(sph_rad.^2 - cyl_rad.^2);
cyl_len = cyl_len-cyl_shift;

sn = 20;
[sx sy sz] = sphere(sn);
ssz = size(sx);

size(sx)
size(sph_rad)
sx = sx * sph_rad;
sy = sy * sph_rad;
sz = sz * sph_rad;

[cx cy cz] = cylinder(cyl_rad*[1 1],sn);
csz = size(cx);
cz = cz*cyl_len + cyl_shift;

ctri = delaunay(cx(2,1:end-1),cy(2,1:end-1));


e = eye(3);
v3 = vec';
v3 = v3/sqrt(v3'*v3);

[m i] = min(sum((repmat(v3,1,3)-e).^2,1));
e23 = e(:,1:3 ~= i);
v1 = cross(v3,e23(:,1));
v1 = v1/sqrt(v1'*v1);

v2 = cross(v1,v3);
v2 = v2/sqrt(v2'*v2);

v = [v1 v2 v3];

C = [cx(:) cy(:) cz(:)]';
S = [sx(:) sy(:) sz(:)]';

C = v*C;
S = v*S;

cx = reshape(C(1,:),csz);
cy = reshape(C(2,:),csz);
cz = reshape(C(3,:),csz);

sx = reshape(S(1,:),ssz);
sy = reshape(S(2,:),ssz);
sz = reshape(S(3,:),ssz);

cx = cx + pos(1);
cy = cy + pos(2);
cz = cz + pos(3);

sx = sx + pos(1);
sy = sy + pos(2);
sz = sz + pos(3);

np = get(gca,'nextplot');
set(gca,'nextplot','add');

prop = {'facecolor','r','linestyle','none',...
    'facelighting','phong',...
	'AmbientStrength',0.3,...
	'DiffuseStrength',0.8,...
	'SpecularStrength',0,...
	'SpecularExponent',10,...
	'SpecularColorReflectance',1,...
    varargin{:}};

surf(sx,sy,sz,prop{:});
surf(cx,cy,cz,prop{:});
trisurf(ctri,cx(2,:),cy(2,:),cz(2,:),prop{:});

set(gca,'nextplot',np);

