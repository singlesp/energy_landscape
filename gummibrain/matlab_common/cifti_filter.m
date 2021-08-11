function newdata = cifti_filter(data,TR_sec,hp_hz,lp_hz,tmpprefix)

if(~exist('tmpprefix','var') || isempty(tmpprefix))
    tmpprefix='';
end

tpDim=2;

outprefix=sprintf('%s_%s',tempname,tmpprefix);


datamean=mean(data,tpDim);
data=bsxfun(@minus,data,datamean);

hp_tr= (1/hp_hz)/(2*TR_sec);
lp_tr= (1/lp_hz)/(2*TR_sec);

BOdimX=size(data,1);  BOdimZnew=ceil(BOdimX/100);  BOdimT=size(data,tpDim);
save_avw(reshape([data ; zeros(100*BOdimZnew-BOdimX,BOdimT)],10,10,BOdimZnew,BOdimT),[outprefix '_fakeNIFTI'],'f',[1 1 1 TR_sec]);

if(lp_hz>0)
    cmd_str=sprintf(['fslmaths ' outprefix '_fakeNIFTI -bptf %f %f ' outprefix '_fakeNIFTI'],hp_tr,lp_tr);
else
    cmd_str=sprintf(['fslmaths ' outprefix '_fakeNIFTI -bptf %f %f ' outprefix '_fakeNIFTI'],hp_tr,-1);
end

system(cmd_str);

grot=reshape(read_avw([outprefix '_fakeNIFTI']),100*BOdimZnew,BOdimT);  

newdata=grot(1:BOdimX,:);  clear grot;

unix(['rm ' outprefix '_fakeNIFTI.nii.gz']);    

newdata=bsxfun(@plus,newdata,datamean);
