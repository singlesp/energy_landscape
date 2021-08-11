function A = mesh_surface_area(faces,verts)

L = reshape(edgelengths(faces,verts),[],3);
s = sum(L,2)/2;
A = sqrt(s.*prod(bsxfun(@minus,s,L),2));
A=sum(A);
%avgvtxarea = sum(A)/size(verts,1);
