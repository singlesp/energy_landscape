function SizeWindow(h, sz)
if(nargin < 1 || isempty(h))
    h = gcf;
end

if(nargin < 2 || isempty(sz))
    wsz = get(0,'ScreenSize');
    sz = wsz([3 4]);
end

pos = get(h,'position');
pos([3 4]) = sz;
set(h,'position',pos);