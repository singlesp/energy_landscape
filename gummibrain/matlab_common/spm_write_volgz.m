function spm_write_vol_gz(Vhdr,vol,fname)

if(exist('fname','var') && ~isempty(fname))
    Vhdr.fname=fname;
else
    fname=Vhdr.fname;
end

[outdir,fnametmp,ext]=fileparts(fname);
if(strcmpi(ext,'.gz'))
    tmp=tempdir;
    Vtmp=rmfield(Vhdr,'pinfo');
    Vtmp.fname=fullfile(tmp,fnametmp);
    spm_write_vol(Vtmp,vol);
    gzip(Vtmp.fname,outdir);
    delete(Vtmp.fname);
else
    spm_write_vol(Vhdr,vol);
end
