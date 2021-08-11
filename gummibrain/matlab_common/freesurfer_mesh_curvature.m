function [curv gausscurv] = freesurfer_mesh_curvature(verts,faces,smoothiter)

if(nargin < 3)
    smoothiter = 0;
end

[tmpd,tmpf,~] = fileparts(tempname);
srffile = [tmpd '/lh.' tmpf];
write_surf(srffile,verts,faces);

if(smoothiter > 0)
    cmd = sprintf('mris_curvature -w -a %d %s',smoothiter,srffile);
else
    cmd = sprintf('mris_curvature -w %s',srffile);
end

[status,result] = system(cmd);

hfile = [srffile '.H'];
kfile = [srffile '.K'];
hcurv = fast_read_curv(hfile);
kcurv = fast_read_curv(kfile);

curv = hcurv;
if(nargout > 1)
    gausscurv = kcurv;
end

delete(srffile,hfile,kfile);