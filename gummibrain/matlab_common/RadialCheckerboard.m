function pict = RadialCheckerboard(num_wedges,num_rings,tilt_angle,use_log)
if(nargin < 1)
    num_wedges = 8;
end

if(nargin < 2)
    num_rings = 8;
end

if(nargin < 3 || isempty(tilt_angle))
    tilt_angle = 0;
end

tilt_angle = tilt_angle+180;

if(nargin < 4)
    use_log = true;
end

radius = 200;
res = 1000;

[x y] = meshgrid(linspace(-radius,radius,res));
r = sqrt(x.^2+y.^2);
if(use_log)
    r = log10(r);
    r = max(r,1);
end

th = mod(atan2(y,x)+tilt_angle*pi/180, 2*pi);
max_th = max(th(:));
max_r = max(r(:));

pict_th = mod(round((th / (2*pi))*num_wedges),2) == 0;
pict_r = mod(fix((r / max_r)*num_rings),2) == 0;
pict = (pict_th & pict_r) | (~pict_th & ~pict_r);
pict = pict*1;
pict(r >= radius) = .5;

%pict = imresize(pict,.25,'bicubic');
