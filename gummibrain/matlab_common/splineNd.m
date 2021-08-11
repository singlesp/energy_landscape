function [zz r] = splineNd(x,z,xx,r)

if(nargin < 4 || isempty(r))
    r = [];
end

xrange = [min(x,[],1); max(x,[],1)];

x0 = bsxfun(@minus,x,xrange(1,:));
x0 = bsxfun(@rdivide,x0,xrange(2,:)-xrange(1));

xx0 = bsxfun(@minus,xx,xrange(1,:));
xx0 = bsxfun(@rdivide,xx0,xrange(2,:)-xrange(1));
% 
%             x0 = (x-min(x))/(max(x)-min(x));
%             y0 = (y-min(y))/(max(y)-min(y));
%             xg0 = (xg-min(xg(:)))/(max(xg(:))-min(xg(:)));
%             yg0 = (yg-min(yg(:)))/(max(yg(:))-min(yg(:)));


if(isempty(r))
    [sp, r] = tpaps(x0.',z.');
else
    sp = tpaps(x0.',z.',r);
end
zz = fnval(sp,xx0.').';
