function [tri vert vals norms] = GetTriSurf(arg1)

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

tri = get(hsurf,'Faces');
vert = get(hsurf,'Vertices');

if(nargout > 2)
    vals = get(hsurf,'FaceVertexCData');
end

if(nargout > 3)
    norms = get(hsurf,'VertexNormals');
    norms = norms ./ repmat(sqrt(sum(norms.^2,2)),1,3);
end