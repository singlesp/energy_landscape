function [rgb_new ratio] = rgblum2rgb(rgb, luminance, displaycalib)

v = .7;
L0 = displaycalib.L0;
L1 = rgb2lum(rgb,displaycalib);
Lv = rgb2lum(v*rgb, displaycalib);

C = [L0 Lv L1];
dispcalC = ComputeCalibration(C,C,C,[0 v 1]);

vnew = ((luminance-L0)/(L1-L0))^(1/dispcalC.bR);
rgb_new = rgb*vnew;

if(nargout > 1)
    ratio = vnew;
end