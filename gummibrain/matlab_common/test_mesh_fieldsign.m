clearb all;
close all;
clc;

if(isunix)
    lapacklib='-lmwlapack';
else
    lapacklib = fullfile(matlabroot,'extern', 'lib', computer('arch'), 'microsoft', 'libmwlapack.lib');
end
mex('-v','-largeArrayDims','mesh_fieldsign_double.c',lapacklib);

%%

wcfile = 'D:\Data\rivalry_eegmri_prism\GROUP\MRI\nii\GROUP_nga._wedge_cos.rh.ico7.w';
wsfile = 'D:\Data\rivalry_eegmri_prism\GROUP\MRI\nii\GROUP_nga._wedge_sin.rh.ico7.w';
rcfile = 'D:\Data\rivalry_eegmri_prism\GROUP\MRI\nii\GROUP_nga._ring_cos.rh.ico7.w';
rsfile = 'D:\Data\rivalry_eegmri_prism\GROUP\MRI\nii\GROUP_nga._ring_sin.rh.ico7.w';

surffile = 'D:\Data\freesurfer\rh.sphere.ico7.reg';

[verts, faces] = freesurfer_read_surf_kj(surffile);
conn = vertex_neighbours(struct('vertices',verts,'faces',faces));
numverts = size(verts,1);

wc = fast_read_wfile_kj(wcfile,numverts);
ws = fast_read_wfile_kj(wsfile,numverts);
rc = fast_read_wfile_kj(rcfile,numverts);
rs = fast_read_wfile_kj(rsfile,numverts);

Veccen = complex(rc,rs);
Vpolar = complex(wc,ws);
pos = verts;
Ne = conn;
eccen_offset = 0;
polar_offset = 0;
smoothiter = 30;

%%
Vfs = mesh_fieldsign(Veccen,Vpolar,pos,Ne,smoothiter,eccen_offset,polar_offset);
%%


figure;
trisurf(faces,verts(:,1),verts(:,2),verts(:,3),sign(Vfs),'linestyle','none');
shading interp;
axis vis3d equal;
