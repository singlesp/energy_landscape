function fwhm = mesh_iter2fwhm(faces,verts,niters)

L = reshape(edgelengths(faces,verts),[],3);
s = sum(L,2)/2;
A = sqrt(s.*prod(bsxfun(@minus,s,L),2));
avgvtxarea = sum(A)/size(verts,1);

%from freesurfer mri_surf2surf.c
%1.14 is a fudge factor based on empirical fit of nearest neighbor
gstd = sqrt((7*avgvtxarea)*(niters - 0.5)/(1.14*4*pi));
fwhm = gstd*sqrt(log(256.0));
