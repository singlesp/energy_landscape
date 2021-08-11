function [verts faces norms] = SimplifyBVSurf(bv, keepratio)

loc = bv.VertexCoordinate;
tri = bv.TriangleVertex;

if(keepratio > 0)
    [loc1, tri1] = meshcheckrepair(loc, tri);
    [loc_simp, tri_simp] = meshresample(loc1,tri1,keepratio);
    verts = bv2tal(loc_simp)';
    faces = tri_simp;

    FV = struct;
    FV.vertices = verts;
    FV.faces = faces;

    norms = patchnormals(FV);
else
    verts = bv2tal(loc)';
    faces = tri;
    norms = bv.VertexNormal;
end