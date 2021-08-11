function imgnew = scaleimage(img,scalefactor,bgcolor,method,kernsize)
%quick version of imresize that does not require image toolbox or opengl

sz=size(img);
sznew=round(sz(1:2)*scalefactor);

[x,y] = meshgrid(linspace(0,1,sz(2)),linspace(0,1,sz(1)));
[x2,y2] = meshgrid(linspace(0,1,sznew(2)),linspace(0,1,sznew(1)));

if(kernsize>0)
    %gaussian is an fsfast function... replace it with something!
    kern=gaussian(linspace(-1,1,kernsize),0,1);
    kern=kern'*kern;
    kern=kern/sum(kern(:));
    kernpad=2*size(kern,1);
    imgpad=padimage(img,kernpad,bgcolor);
end

imgnew=zeros([sznew 3]);

for i = 1:3
    if(kernsize>0)
        imgf=conv2(double(imgpad(:,:,i)),kern,'same');
        imgf=imgf(kernpad+(1:sz(1)),kernpad+(1:sz(2)));
    else
        imgf=double(img(:,:,i));
    end
    imgnew(:,:,i)=interp2(x,y,imgf,x2,y2,method);
end

imgnew=uint8(min(255,max(0,round(imgnew))));
