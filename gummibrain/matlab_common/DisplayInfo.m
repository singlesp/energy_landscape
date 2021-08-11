%viewprop = DisplayInfo(window, viewdist_cm, [screensize_cm])
%
%return a bunch of useful information about the display:
%resolution, physical width/height, ppd, fov etc...
%
%NOTE: the EDID monitor information seems to provide inaccurate physical
% size much of the time.  Would be better off measuring and inputing these
% manually in screensize_cm as [width_cm height_cm]
function viewprop = DisplayInfo(window, viewdist_cm, screensize_cm)

[wpx hpx] = Screen('WindowSize',window);

if(nargin < 3 || isempty(screensize_cm))
    [wmm hmm] = Screen('DisplaySize',window); 
    wcm = wmm/10;
    hcm = hmm/10;
else
    wcm = screensize_cm(1);
    hcm = screensize_cm(2);
end

ax = 2*atan2(wcm/2,viewdist_cm)*180/pi;
ay = 2*atan2(hcm/2,viewdist_cm)*180/pi;

viewprop = struct;
viewprop.width_cm = wcm;
viewprop.height_cm = hcm;
viewprop.distance_cm = viewdist_cm;
viewprop.size_cm = [wcm hcm];

viewprop.width = wpx;
viewprop.height = hpx;
viewprop.size = [wpx hpx];
viewprop.center = round(viewprop.size/2);
viewprop.ppd = wpx/ax;
viewprop.ppcm = wpx/wcm;
viewprop.ppdxy = [wpx/ax hpx/ay];
viewprop.ppcmxy = [wpx/wcm hpx/hcm];
viewprop.fovx = ax;
viewprop.fovy = ay;
