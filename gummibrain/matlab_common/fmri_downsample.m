function fmri_downsample(niifile,outfile,dsfactor)

M=load_nifti(niifile);
volshape=size(M.vol);

S=reshape(M.vol,[],size(M.vol,4)).';

Smean=mean(S,1);
S=bsxfun(@minus,S,Smean);

dsfactor_orig=dsfactor;
ds_p=1;
ds_q=dsfactor;
if(dsfactor ~= floor(dsfactor))
    [ds_q,ds_p]=rat(dsfactor, .01);
    dsfactor=ds_q/ds_p;
end

fprintf('\ninput dsfactor = %g\nq/p = %d/%d\nnew dsfactor = %g\n',...
    dsfactor_orig,ds_q,ds_p,dsfactor);
    
if(ds_p > 15)
    %do we need this limit?  at some point it blows up right?
   %error('upsampling capped at 15x.  comment out this error to continue.');
end

Snew=resample(S,ds_p,ds_q);
Snew=bsxfun(@plus,Snew,Smean);
Snew=reshape(Snew.',[volshape(1:3),size(Snew,1)]);

Mnew=M;
Mnew.vol=Snew;
Mnew.dim(5)=size(Snew,2);
Mnew.pixdim(5)=M.pixdim(5)*dsfactor;

save_nifti(Mnew,outfile);

%saves with pixdim4 in msec, but some programs assume sec
%so use fslchfiletype to just canonicalize


[d,f,ext]=fileparts(outfile);
if(strcmpi(ext,'.nii.gz'))
    system(['fslchfiletype NIFTI_GZ ' fullfile(d,f)]);
elseif(strcmpi(ext,'.nii'))
    system(['fslchfiletype NIFTI ' fullfile(d,f)]);
end
system(['fslmaths ' fullfile(d,f) ' ' fullfile(d,f)]);

if(strcmpi(ext,'.nii.gz'))
    system(['rm -f ' fullfile(d,[f '.nii'])]);
elseif(strcmpi(ext,'.nii'))
    system(['rm -f ' fullfile(d,[f '.nii.gz'])]);
end
