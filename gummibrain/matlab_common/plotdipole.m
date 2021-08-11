function [hcylinder hsphere] = plotdipole(pos,vec,cvalues,varargin)

show_sphere = true;
show_cyl = true;

cylres = 30;
scalefactor = 1;
cylfactor = .5;
cylaspect = 3;

if(nargin > 3 && isnumeric(varargin{1}))
    scalefactor = varargin{1};
    if(nargin > 3)
        varargin = {varargin{2:end}};
    else
        varargin = {};
    end
end

varargin = {'edgealpha',0,varargin{:}};
if(nargin < 3 || isempty(cvalues))
    cvalues = [];
end
cvalues = cvalues(:); %make sure it's a column vector


N = size(pos,1);

if(numel(cvalues) == 1)
    cvalues = cvalues*ones(N,1);
end

D = vec;
Dlen = sqrt(sum(D.^2,2));
vecsize = Dlen;

horz_idx = D(:,3) == 0;

Daxis = cross(D,[D(:,1:2) zeros(size(D,1),1)],2);
Daxis_horz = cross(D,[D(:,1:2) -ones(size(D,1),1)],2);

Daxis(horz_idx,:) = Daxis_horz(horz_idx,:); %need other crossprod reference for horz lines
Daxis_len = sqrt(sum(Daxis.^2,2));

%Daxis(Daxis_len == 0,:) = repmat([1 0 0],sum(Daxis_len == 0), 1);
Daxis_len(Daxis_len == 0) = 1;
Daxis = Daxis./repmat(Daxis_len,[1 3]);


Dtheta = sign(D(:,3)).*acos(D(:,3)./Dlen);
Dtheta(horz_idx) = acos(0); %horz lines are always pi;
%%
hsph = [];
hcyl = [];


if(show_sphere)
    hsph = zeros(size(pos,1),1);
end

if(show_cyl)
    hcyl = zeros(size(pos,1)-1,1);
end

%if(show_sphere)
%    hsph = plotsphere(pos,vecsize*scalefactor,cvalues,'resolution',cylres,varargin{:});
%end

[Xs Ys Zs] = sphere(cylres);
[X Y Z] = cylinder(1,cylres);

%X = X-mean(X(:));
np = get(gca,'nextplot');
set(gca,'nextplot','add');

for i = 1:size(pos,1)
    xi = cylfactor*scalefactor*vecsize(i).*X;
    yi = cylfactor*scalefactor*vecsize(i).*Y;
    zi = cylaspect*scalefactor*vecsize(i).*Z;
    
    xsi = scalefactor*vecsize(i)*Xs;
    ysi = scalefactor*vecsize(i)*Ys;
    zsi = scalefactor*vecsize(i)*Zs;

    cp = [xi(2,:); yi(2,:); zi(2,:)]';
    ctri = delaunay(cp(1:end-1,1), cp(1:end-1,2));

    cyltri = repmat([1:numel(xi)-2]',[1 3]);
    cyltri(:,2) = cyltri(:,1) + 1;
    cyltri(:,3) = cyltri(:,1) + 2;

    cyltri(2:2:end, [1 3 2]) = cyltri(2:2:end,[1 2 3]);
    cyltri = [cyltri; numel(xi)+1+ctri(:,[1 3 2])];
    
    %xsi = spherescale*width(i)*Xs;
    %ysi = spherescale*width(i)*Ys;
    %zsi = spherescale*width(i)*Zs;

    [ptmp R] = RodrigRotate(Daxis(i,:),Dtheta(i),[0 0 1]);
    
    p = [R*[xi(:) yi(:) zi(:)]']';
    p = p + repmat(pos(i,:),size(p,1),1);
    xi = reshape(p(:,1),2,size(p,1)/2);
    yi = reshape(p(:,2),2,size(p,1)/2);
    zi = reshape(p(:,3),2,size(p,1)/2);
    
    p = [R*[xsi(:) ysi(:) zsi(:)]']';
    p = p + repmat(pos(i,:),size(p,1),1);
    xsi = reshape(p(:,1),size(xsi));
    ysi = reshape(p(:,2),size(ysi));
    zsi = reshape(p(:,3),size(zsi));
    
    tp = [xi(:) yi(:) zi(:)];
    tp = [tp; tp(2:2:end,:)];

    cv = [];
    if(~isempty(cvalues))
        cv = cvalues(mod(i-1,numel(cvalues))+1);
    end
    
    if(show_cyl)
        hcyl(i) = trisurf(cyltri,tp(:,1),tp(:,2),tp(:,3),varargin{:}); 
        if(~isempty(cv))
            set(hcyl(i),'CData',cv);
        end
    end

    if(show_sphere)
        hsph(i) = surf(xsi,ysi,zsi,varargin{:});
        if(~isempty(cv))
            set(hsph(i),'CData',repmat(cv,size(Xs)));
        end
    end
end

set(gca,'nextplot',np);

if(nargout > 0)
    hcylinder = hcyl;
end

if(nargout > 1)
    hsphere = hsph;
end