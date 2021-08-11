function [surf_val valid_idx hsurf] = GetVertexColors(arg1)
%returns the vertex values of a surface/patch object.  This can be either 
%numerical or RGB, depending on the surface.  Since some surfaces contain
%vertices that are not part of valid faces/triangles, this function can
%also return a list of which vertices actually belong to faces
%
%surf_val = GetVertexColors
%surf_val = GetVertexColors(ax)
%surf_val = GetVertexColors(hsurf)
%[surf_val valid_idx] = GetVertexColors(ax)
%
%
%ax = axis in which to find the surface (default = gca)
%hsurf = surface object to get vertices for
%
%surf_val = Nx1 or Nx3 matrix with the values from each vertex (valid or not)
%valid_idx = Mx1 list containing a subset of values 1:N, indicating which
%   of the N vertices belong to valid faces.

if(nargin < 1)
    arg1 = gca;
end

try 
    switch lower(get(arg1,'type'))
        case 'axes'
            ax = arg1;
            ch = get(ax,'children');
            hsurf = ch(strcmpi(get(ch,'type'),'patch'));
        case 'patch'
            hsurf = arg1;
        otherwise
            error('input is not an axes or patch object');
    end
catch    
    error('input is not an axes or patch object');
end

surf_val = get(hsurf,'FaceVertexCData');

if(nargout > 1)
    surf_tri = get(hsurf,'Faces');
    valid_idx = unique(surf_tri(:));
end