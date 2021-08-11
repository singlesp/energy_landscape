function [handles valididx] = CreateTopoSubplots(locs,locmask,masktype,subwidth,subheight)

handles = [];

%%% read .locs file if filename given
if(ischar(locs))
    L = readlocs(locs);
    locs = L;
end

%%% if not specified, assume all electrodes should be shown
if(nargin < 2 || isempty(locmask))
    locsmask = 1;
end

if(nargin < 3 || isempty(masktype))
    masktype = true;
end


if(iscell(locmask))
    locmask = StringMask({locs.labels},locmask,masktype);
end

if(numel(locmask) == 1)
    locmask = locmask*ones(size(locs));
end

if(nargin < 4)
    subwidth = 1/20;
    subheight = 1/15;
end

th = [locs.theta]'*pi/180;
r = [locs.radius]';
r = .5*r/max(r);

y = .5 + r.*cos(th) - subwidth/2;
x = .5 + r.*sin(th) - subheight/2;
pos = [x y subwidth*ones(size(x)) subheight*ones(size(y))];
pos(~locmask,1:2) = -1;
pos(~locmask,3:4) = 0;

handles = zeros(numel(locs),1);

for i = 1:numel(locs)
    if(~locmask(i))
        continue;
    end

    h = subplot('position',[-1 -1 subwidth subheight]);
    %set(h,'outerposition',pos(i,:));
    set(h,'position',pos(i,:));

    handles(i) = h;


end

if(nargout > 1)
    valididx = find(locmask);
end

