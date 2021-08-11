function V = rotatevolume(V0,az_el_tilt,interptype,scalefactor)
%V = rotatevolume(V0,az_el_tilt,interptype,scalefactor)
%V = rotatevolume(V0,az_el_tilt,scalefactor)


viewazel=[0 0 0];
viewazel(1:numel(az_el_tilt))=az_el_tilt;

default_interptype='cubic';

if(~exist('interptype','var') || isempty(interptype))
    interptype=default_interptype;
end

if(~exist('scalefactor','var') || isempty(scalefactor))
    scalefactor=1;
end

if(ischar(interptype))
    if(ismember(lower(interptype),{'nearest','linear','cubic','spline'}))
        interptype=lower(interptype);
    elseif(strcmpi(interptype,'trilinear'))
        interptype='linear';
    else
        error('invalid interp type: %s',interptype);
    end
elseif(~isempty(interptype) && isnumeric(interptype))
    scalefactor=interptype;
    interptype=default_interptype;
end

do_square=false;
if(do_square)
    sz0=size(V0);
    sz=max(sz0)*[1 1 1];
    V=zeros(sz(1),sz(2),sz(3));
    V(1:sz0(1),1:sz0(2),1:sz0(3))=V0;
    V=circshift(V,round((sz-sz0)/2));
else
    V=V0;
    sz=size(V);
end


sz2=round(sz*scalefactor);

mx=linspace(1,sz(1),sz2(1));
my=linspace(1,sz(2),sz2(2));
mz=linspace(1,sz(3),sz2(3));
[gy,gx,gz]=meshgrid(my,mx,mz);
xyz=[gx(:) gy(:) gz(:)];


%M=RotationMatrix_azel(180-viewazel(1),viewazel(2),viewazel(3));
%M=RotationMatrix_azel(viewazel(1),viewazel(2),viewazel(3));
M=RotationMatrix(-viewazel(2)*pi/180,-viewazel(3)*pi/180,-viewazel(1)*pi/180);
%M=round(M*1e8)/1e8
Mshift=eye(4);
Mshift(1:3,4)=sz/2;

M=Mshift*M*inv(Mshift);

xyz2=affine_transform(M,xyz);

try
    %use GPU for rotation interp if available (much faster)
    isgpu=gpuDeviceCount;
    if(isgpu>0)
        V=reshape(gather(interpn(gpuArray(double(V)),xyz2(:,1),xyz2(:,2),xyz2(:,3),interptype)),sz2);
    end
catch
    isgpu=0;
end

if(isgpu==0)
    %V=reshape(interp3(double(V),xyz2(:,2),xyz2(:,1),xyz2(:,3),interptype),round(sz*scalefactor));
    G=griddedInterpolant(double(V),interptype); %faster than interp3
    V=reshape(G(xyz2(:,1:3)),sz2); %note [x y z] instead of interp3(y,x,z)
end

V(~isfinite(V))=0;
V=cast(V,'like',V0);