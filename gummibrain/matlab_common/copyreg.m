function newfile = copyreg(fromfile,tofile,prefix,fileext)
if(~exist('prefix','var') || isempty(prefix))
    prefix = 'copyreg_';
end
if(~exist('fileext','var'))
    fileext = '.nii';
end


if(~iscell(tofile))
    tofile = {tofile};
end

newfile = {};
mat = spm_get_space(fromfile);
for f = 1:numel(tofile)
    [fd,fn,fext] = fileparts(tofile{f});
    if(~isempty(fileext))
        fext = fileext;
    end
    newfile{f} = sprintf('%s/%s%s%s',fd,prefix,fn,fext);
    Vhdr = spm_vol(tofile{f});
    V = spm_read_vols(Vhdr);
    Vhdr.fname = newfile{f};
    spm_write_vol(Vhdr,V);
    spm_get_space(newfile{f},mat);
end
