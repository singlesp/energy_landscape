function newfile = mri_file(fname)

do_outchar = ischar(fname);
if(do_outchar)
    fname = cellstr(fname);
end

fname = regexprep(fname,'\.hdr$','.img');

isimg = regexpimatch(fname,'\.img$');
isnii = regexpimatch(fname,'\.nii$');
isblank = ~(isimg | isnii);


newnii = fname;
newnii(isblank) = strcat(fname(isblank),'.nii');
nii_exist = cellfun(@(f)(exist(f,'file')>0),newnii);

newimg = fname;
newimg(isblank) = strcat(fname(isblank),'.img');
img_exist = cellfun(@(f)(exist(f,'file')>0),newimg);

newfile = fname;
newfile(isblank & img_exist) = newimg(isblank & img_exist);
newfile(isblank & nii_exist) = newnii(isblank & nii_exist);

if(do_outchar)
    newfile = char(newfile);
end