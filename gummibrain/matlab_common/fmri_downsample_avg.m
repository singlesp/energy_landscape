function fmri_downsample_avg(niifile,outfile,dsfactor)

M=load_nifti(niifile);
%M=load_nifti(niifile,true);
%[M.vol,dims,scales,bpp,endian]=read_avw(niifile);

volshape=size(M.vol);

S=reshape(M.vol,[],size(M.vol,4)).';

Smean=mean(S,1);
S=bsxfun(@minus,S,Smean);

%Snew=resample(S,1,dsfactor);
skipidx=1:dsfactor:(size(S,1)-dsfactor+1);
Snew=zeros(numel(skipidx),size(S,2));
for i = 1:dsfactor
    Snew=Snew+S(skipidx+i-1,:);
end
Snew=Snew/dsfactor;
Snew=bsxfun(@plus,Snew,Smean);
Snew=reshape(Snew.',[volshape(1:3),size(Snew,1)]);

Mnew=M;
Mnew.vol=Snew;
Mnew.dim(5)=size(Snew,2);
Mnew.pixdim(5)=M.pixdim(5)*dsfactor;

save_nifti(Mnew,outfile);
%scales(4)=scales(4)*dsfactor;
%save_avw(Mnew.vol,outfile,'f',scales);
%return
 
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
