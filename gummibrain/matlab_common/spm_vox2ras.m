function [m_vox2ras m_cras] = spm_vox2ras(mri)
% [m_vox2ras m_cras] = spm_vox2ras(mri)
%
% Return freesurf style vox2ras affine transform and c_ras transform
%   so m_cras*(freesurf verts) is in voxel space for coregistration or whatever
%
if(ischar(mri) && exist(mri,'file'))
    V = spm_vol(mri);
elseif(isstruct(mri))
    V = mri;
end

%m_vox2ras = [V.private.hdr.srow_x; V.private.hdr.srow_y; V.private.hdr.srow_z; 0 0 0 1];
m_vox2ras = V.mat + [zeros(3) sum(V.mat(1:3,1:3),2); 0 0 0 0];
m_cras = spm_matrix(affine_transform(m_vox2ras,V.dim/2));
