function convert_stl_to_freesurfer(stlfilename,outfilename,reffilename)
[vnew,fnew]=stlread(stlfilename);
[v,f,tag]=freesurfer_read_surf_kj(reffilename);
freesurfer_write_surf_kj(outfilename,vnew,fnew,tag);
