function img=surf2volz_kj(node,face,xi,yi,zi)
%
% img=surf2volz(node,face,xi,yi,zi)
%
% convert a triangular surface to a shell of voxels in a 3D image
% along the z-axis
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% input:
%	 node: node list of the triangular surface, 3 columns for x/y/z
%	 face: triangle node indices, each row is a triangle
%	 xi,yi,zi: x/y/z grid for the resulting volume
%
% output:
%	 img: a volumetric binary image at position of ndgrid(xi,yi,zi)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

ne=size(face,1);
%img=zeros(length(xi),length(yi),length(zi),'uint8');

dx0=min(abs(diff(xi)));
dx=dx0/2;
dy0=min(abs(diff(yi)));
dy=dy0/2;
dz0=min(abs(diff(zi)));

dl=sqrt(dx*dx+dy*dy);

minz=min(node(:,3));
maxz=max(node(:,3));
iz=hist([minz,maxz],zi);
hz=find(iz);
iz=hz(1):min(length(zi),hz(end)+1);

lenman = abs(dx)+abs(dy);
z=zi(iz);

imgc = cell(length(z),1);
zz = zeros(length(z),1);


edges=[ face(:,[1,2]); face(:,[1,3]); face(:,[2,3]) ];

for i=1:length(z)
    
    plane = z(i);
    [bcutpos,bcutedges]=qmeshcutz_kj(size(face),node,plane,edges);

    if(isempty(bcutpos))
        continue;
    end

    e0 = bcutpos(bcutedges(:,1),1:2);
    e1 = bcutpos(bcutedges(:,2),1:2);
    dd=e1-e0;
    len = ceil(sum(abs(dd),2)/lenman)+1;
    dd=dd./[len len];
    posx = floor(bsxfun(@plus,e0(:,1),dd(:,1)*(0:len)-xi(1))/dx0)';
    posy = floor(bsxfun(@plus,e0(:,2),dd(:,2)*(0:len)-yi(1))/dy0)';
    vcheck = ~(posx>length(xi) | posy>length(yi) | posx<=0 | posy<=0);
    pos = [posx(vcheck(:)) posy(vcheck(:))]; 

    if(isempty(pos))
        continue;
    end

    zz(i)=floor(((z(i)-zi(1)))/dz0);
    imgc{i} = zeros(length(xi),length(yi));
    posidx = sub2ind(size(imgc{i}),pos(:,1),pos(:,2));
    imgc{i}(posidx) = 1;
    
    
end

img = zeros(length(xi),length(yi),length(zi),'uint8');
img(:,:,zz(zz>0)) = cat(3,imgc{zz>0});
