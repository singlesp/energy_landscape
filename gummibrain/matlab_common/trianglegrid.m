function [tri tx ty] = trianglegrid(x,y)

if(~isvector(x))
    x = x(:);
end
if(~isvector(y))
    y = y(:);
end

nx = numel(x);
ny = numel(y);
n = nx*ny;

%tri = zeros(2*(numel(x)-1)*(numel(y)-1),3);

ta = [(1:nx-1)' (2:nx)' (nx+1:2*nx-1)'];
tb = [(nx+1:2*nx-1)' (2:nx)' (nx+2:2*nx)'];
tri = [];
for j = 1:ny-1
    tri = [tri; ta + (j-1)*nx; tb + (j-1)*nx];
end

if(nargout > 1)
    [ty tx] = meshgrid(x,y);
    tx = tx(:);
    ty = ty(:);
end