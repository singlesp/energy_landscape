function img = fullcolors(width, height, borderheight)

sz = [width height];
bordersz = [width borderheight];

[x y] = meshgrid(linspace(0,1,sz(1)),linspace(0,1,sz(2)/2));
hsvimg = zeros([sz([2 1]) 3]);

hsvimg(:,:,1) = [x; x];
hsvimg(:,:,2) = [y; ones(size(y))];
hsvimg(:,:,3) = [ones(size(x)); 1-y];

img = [ones([bordersz([2 1]) 3]); 
    hsv2rgb(hsvimg); 
    zeros([bordersz([2 1]) 3])];