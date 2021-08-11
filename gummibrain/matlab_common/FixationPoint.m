function img = FixationPoint(sz, thickness, tilt_angle)

if(nargin < 2 || isempty(thickness))
    thickness = 1;
end

thickness = thickness-1;
img = ones(sz,sz);
img(:,round(sz/2 - thickness/2):(round(sz/2 + thickness/2))) = 0;
img(round(sz/2 - thickness/2):(round(sz/2 + thickness/2)),:) = 0;