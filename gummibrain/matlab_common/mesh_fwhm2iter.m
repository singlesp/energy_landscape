function niters = mesh_fwhm2iter(faces,verts,fwhm)

L = reshape(edgelengths(faces,verts),[],3);
s = sum(L,2)/2;
A = sqrt(s.*prod(bsxfun(@minus,s,L),2));
avgvtxarea = sum(A)/size(verts,1);

%from freesurfer mri_surf2surf.c
%1.14 is a fudge factor based on empirical fit of nearest neighbor
gstd = fwhm/sqrt(log(256.0));
niters = floor(1.14*(4*pi*(gstd*gstd))/(7*avgvtxarea) + 0.5);
