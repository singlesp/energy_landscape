function img = graycolors(width, height, borderheight)

sz = [width height];
bordersz = [width borderheight];

rgbimg = repmat(linspace(1,0,sz(2))',[1 sz(1) 3]);

img = [ones([bordersz([2 1]) 3]); 
    rgbimg; 
    zeros([bordersz([2 1]) 3])];