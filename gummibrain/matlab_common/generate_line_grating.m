function img = generate_line_grating(theta, spatial_freq, linewidth, imgspacing, imgwidth)
%theta = angle in radians
%spatial_freq = cycles/image coordinate
%imgspacing = point spacing (spatial sampling) in image space (eg, 0.1)
%imgwidth = width of image

[x y] = meshgrid(0:imgspacing:imgwidth);

u = sin(theta);
v = cos(theta);
% u = 1;
% v = 0;
% 
% r = sqrt(u^2+v^2); %normalize
% u = u/r;
% v = v/r;



f = 1/(spatial_freq*imgspacing);
img = cos(2*pi*f*(u*x + v*y));

if(linewidth > 0)
    w = linewidth*imgspacing/2;
    cutoff = cos(2*pi*f*w);
    img = img > cutoff;
end