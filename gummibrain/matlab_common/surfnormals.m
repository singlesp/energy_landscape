% surf_norms = surfnormals(verts,tri)
% verts = Nx3
% tri = Mx3
function surf_norms = surfnormals(verts, tri)

FV = struct;
FV.vertices = verts;
FV.faces = tri;
surf_norms = patchnormals(FV);