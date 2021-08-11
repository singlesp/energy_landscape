function fmri_downsample(niifile,outfile,dsfactor,newTR)

if(~exist('newTR','var'))
    newTR=0;
end 

wbc='/home/range1-raid1/kjamison/workbench_v1.2.3/bin_rh_linux64/wb_command';
BO=ciftiopen(niifile,wbc);
%M=load_nifti(niifile);
%volshape=size(M.vol);

%S=reshape(M.vol,[],size(M.vol,4)).';

S=double(BO.cdata).';

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
%Snew=reshape(Snew.',[volshape(1:3),size(Snew,1)]);

BO.cdata=single(Snew.');

ciftisave(BO,outfile,wbc);
if(newTR ~= 0)
    system(sprintf('%s -cifti-change-mapping %s ROW %s -series %.3f 0 -unit SECOND',wbc,outfile,outfile,newTR));
end
