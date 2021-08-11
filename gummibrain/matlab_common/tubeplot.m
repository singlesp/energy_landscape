function [hcylinder hsphere] = tubeplot(verts,widths,cvalues,args,varargin)


%rotate_sphere = false;

if(nargin < 4 || isempty(args))
    argnames = [];
    argvals = [];
    argres = -1;
    arginterp = -1;
    argsph = -1;
    argcyl = -1;
    argsphscale = -1;
else
    argnames = {args{1:2:end}};
    argvals = {args{2:2:end}};
    argres = StringIndex(argnames,'resolution');
    arginterp = StringIndex(argnames,'interp');
    argsph = StringIndex(argnames,'spheres');
    argcyl = StringIndex(argnames,'cylinders');
    argsphscale = StringIndex(argnames,'spherescale');
end


if(argres > 0)
    cylres = argvals{argres};
else
    cylres = 20;
end

if(arginterp > 0)
    interpfactor = argvals{arginterp};
else
    interpfactor = 1;
end

if(argsph > 0)
    show_sphere = argvals{argsph};
else
    show_sphere = true;
end

if(argcyl > 0)
    show_cyl = argvals{argcyl};
else
    show_cyl = true;
end

if(argsphscale > 0)
    spherescale = argvals{argsphscale};
else
    spherescale = .97;
end

if(nargin < 3 || isempty(cvalues))
    cvalues = [];
end
cvalues = cvalues(:); %make sure it's a column vector


N = size(verts,1);

if(numel(widths) == 1)
    widths = widths*ones(N,1);
end

verts = spline(linspace(0,1,N),verts',linspace(0,1,interpfactor*N))';
widths = spline(linspace(0,1,N),widths',linspace(0,1,interpfactor*N));
cvalues = spline(linspace(0,1,N),cvalues',linspace(0,1,interpfactor*N));

D = diff(verts);
Dlen = sqrt(sum(diff(verts).^2,2));

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
    
    hsph = plotsphere(verts,widths*spherescale,cvalues,'resolution',cylres,varargin{:});
end

%[Xs Ys Zs] = sphere(cylres);
[X Y Z] = cylinder(1,cylres);
%X = X-mean(X(:));
hold on;

if(show_sphere)
    hsph = zeros(size(verts,1),1);
end

if(show_cyl)
    hcyl = zeros(size(verts,1)-1,1);
end

for i = 1:size(verts,1)-1
    xi = diag(widths(i:i+1))*X;
    yi = diag(widths(i:i+1))*Y;
    zi = diag([0 Dlen(i)])*Z;
    
    %xsi = spherescale*width(i)*Xs;
    %ysi = spherescale*width(i)*Ys;
    %zsi = spherescale*width(i)*Zs;

    [ptmp R] = RodrigRotate(Daxis(i,:),Dtheta(i),[0 0 1]);
    
    p = [R*[xi(:) yi(:) zi(:)]']';
    p = p + repmat(verts(i,:),size(p,1),1);
    xi = reshape(p(:,1),2,size(p,1)/2);
    yi = reshape(p(:,2),2,size(p,1)/2);
    zi = reshape(p(:,3),2,size(p,1)/2);
    
    %ps = [R*[xsi(:) ysi(:) zsi(:)]']';
    %ps = ps + repmat(verts(i,:),size(ps,1),1);
    %xsi = reshape(ps(:,1),size(xsi));
    %ysi = reshape(ps(:,2),size(ysi));
    %zsi = reshape(ps(:,3),size(zsi));
    
    if(show_cyl)
        hcyl(i) = surface(xi,yi,zi,varargin{:});
        if(~isempty(cvalues))
            set(hcyl(i),'CData',repmat(cvalues(i:i+1)',1,cylres+1));
        end
    end
    %hsph(i) = surface(xsi,ysi,zsi,'facecolor','r','facealpha',1,'edgealpha',0);
end

if(nargout > 0)
    hcylinder = hcyl;
end

if(nargout > 1)
    hsphere = hsph;
end