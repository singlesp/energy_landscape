function hsurf = SetVertexColors(arg2,vals)
%sets the vertex values of a surface/patch object.  This can be either 
%numerical or RGB, depending on the surface.
%
%SetVertexColors(surf_val)
%SetVertexColors(ax,surf_val)
%SetVertexColors(hsurf, surf_val)
%
%
%surf_val = Nx1 or Nx3 matrix with the values to use for each vertex.
%   Every vertex must have a value, whether it belongs to a valid face or not.
%ax = axis in which to find the surface (default = gca)
%hsurf = surface object to set vertices for

if(nargin < 2)
    vals = arg2;
    arg2 = gca;
end

try 
    switch lower(get(arg2,'type'))
        case 'axes'
            ax = arg2;
            ch = get(ax,'children');
            hsurf = ch(strcmpi(get(ch,'type'),'patch'));
        case 'patch'
            hsurf = arg2;
        otherwise
            error('input is not an axes or patch object');
    end
catch    
    error('input is not an axes or patch object');
end

set(hsurf,'FaceVertexCData',vals);
