function gray2fun=GrayToLRSub(graycount,templatedir)
%%% Given a grayordinate data file (eg: 91282 values), how do we map these to
%%% cortical hemispheres, and to ROIs within those hemispheres

if(~exist('graycount','var') || isempty(graycount))
    graycount=91282;
    surfname='32k';
end
if(~exist('templatedir','var') || isempty(templatedir))
    templatedir=[homedir '/Source/Pipelines/global/templates'];
end


if(graycount==91282)
    
else
    error('Not sure how to process %s grayord',graycount);
end


templatename=sprintf('%d_Greyordinates',graycount);
Gstruct=ciftiopen(sprintf('%s/%s/%s.dscalar.nii',templatedir,templatename,templatename));
Lstruct=gifti(sprintf('%s/%s/L.atlasroi.%s_fs_LR.shape.gii',templatedir,templatename,surfname));
Rstruct=gifti(sprintf('%s/%s/R.atlasroi.%s_fs_LR.shape.gii',templatedir,templatename,surfname));



SCmask=Gstruct.cdata>1;
Lmask=Lstruct.cdata>0;
Rmask=Rstruct.cdata>0;

Ridx=find(Rmask>0);
Lidx=find(Lmask>0);

%Ridx_gray=1:numel(Ridx);
%Lidx_gray=(1:numel(Lidx))+numel(Ridx);
Lidx_gray=1:numel(Lidx);
Ridx_gray=(1:numel(Ridx))+numel(Lidx);
SCidx_gray=find(SCmask>0);


g2left=@(G)(subsasgn(nan(size(Lmask,1),size(G,2)),struct('type','()','subs',{{Lidx,':'}}),G(Lidx_gray,:)));
g2right=@(G)(subsasgn(nan(size(Rmask,1),size(G,2)),struct('type','()','subs',{{Ridx,':'}}),G(Ridx_gray,:)));
g2subcort=@(G)(G(SCidx_gray));

left2g=@(L)(subsasgn(nan(size(Gstruct.cdata,1),size(L,2)),struct('type','()','subs',{{Lidx_gray,':'}}),L(Lidx,:)));
right2g=@(R)(subsasgn(nan(size(Gstruct.cdata,1),size(R,2)),struct('type','()','subs',{{Ridx_gray,':'}}),R(Ridx,:)));
subcort2g=@(S)(subsasgn(nan(size(Gstruct.cdata,1),size(S,2)),struct('type','()','subs',{{SCidx_gray,':'}}),S));


gray2fun=struct('gray2left',g2left,'gray2right',g2right,'gray2subcort',g2subcort,...
    'left2gray',left2g,'right2gray',right2g,'subcort2gray',subcort2g,...
    'graycount',graycount);
