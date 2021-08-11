function Vnew = floodfill3d(V,p,thresh)

if(~exist('thresh','var') || isempty(thresh))
    thresh=mean(V(:));
end

Vbin=V>=thresh;
Vnew=(+Vbin + ~imfill(~Vbin,p,26))==1;
