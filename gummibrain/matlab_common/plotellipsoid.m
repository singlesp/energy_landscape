function [hsurf fv] = plotellipsoid(pos,radius,zdir,cvalues,varargin)

if(nargin < 3 || isempty(zdir))
    zdir = [0 0 1];
end

if(nargin < 4 || isempty(cvalues))
    cvalues = [];
end

if(numel(varargin) > 1 && strcmpi(varargin{1},'resolution'))
    trires = varargin{2};
    varargin = {varargin{3:end}};
else
    trires = .1;
end

if(min(size(pos)) == 1)
	pos = reshape(pos,1,3);
end

[sphvert, sphface] = meshunitsphere(trires);
[sphvert, sphface] = meshcheckrepair(sphvert,sphface);
        
np = get(gca,'nextplot');
set(gca,'nextplot','add');
hsurf = zeros(size(pos,1),1);

e3 = zdir;
e1 = normrows(([0 -1 0; 1 0 0; 0 0 0]*e3')');
e2 = normrows(cross(e1,e3));

fv = struct();
M = zeros(4,4);
for i = 1:size(pos,1)
    nr = min(i,size(radius,1));
    nz = min(i,size(e3,1));
    
    M(:) = 0;
    M(4,4) = 1;
    M(1:3,1:3) = eye(3).*diag(radius(nr,:));
    if(~all(e3(nz,:) == [0 0 1]))
        M(1:3,1:3) = [e1(nz,:); e2(nz,:); e3(nz,:)]'*M(1:3,1:3);
    end
    M(1:3,4) = pos(i,:)';
    v = affine_transform(M,sphvert);
    
    
    fv(i).faces = sphface;
    fv(i).vertices = v;
    
    hsurf(i) = trisurf(sphface,v(:,1),v(:,2),v(:,3),varargin{:});
    
    if(~isempty(cvalues))
        nc = min(i,size(cvalues,1));
        set(hsurf(i),'FaceColor',cvalues(nc,:));
    end
        
end
set(gca,'nextplot',np);

