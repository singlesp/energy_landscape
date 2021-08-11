function varargout = cloud2surf(cloud)

%switch dims 1 and 2 because of how isosurface expects inputs
cloud(:,[1 2 3]) = cloud(:,[2 1 3]);

volsz = max(cloud)+2;

ind = sub2ind(volsz,cloud(:,1),cloud(:,2),cloud(:,3));

vol = zeros(volsz);
vol(ind) = 1;

fv = isosurface(vol,.5);

if(nargout == 1)
    varargout{1} = fv;
elseif(nargout == 2)
    varargout{1} = fv.faces;
    varargout{2} = fv.vertices;
end

