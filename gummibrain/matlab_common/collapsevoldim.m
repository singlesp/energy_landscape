function mipimg=collapsevoldim(X,bgcolor,mipfun,dims)

if(~exist('bgcolor','var') || isempty(bgcolor))
    bgcolor=0;
end
if(~exist('mipfun','var') || isempty(mipfun))
    mipfun='maxabs';
end
if(~exist('dims','var') || isempty(dims))
    dims=1:3;
end
mipimg={};
sz=size(X);

for di = 1:numel(dims)
    d=dims(di);
    Xtmp=reshape(permute(X,[setdiff(1:3,d) d]),[],size(X,d));

    if(iscell(mipfun))
        mipfun1=mipfun{min(di,numel(mipfun))};
    else
        mipfun1=mipfun;
    end
    if(isnumeric(mipfun1))
        mipslice=mipfun1;
        if(mipslice>=0 && mipslice<1)
            mipslice=round(mipslice*size(Xtmp,2));
        end
        mipslice=max(1,min(size(Xtmp,2),mipslice));
        mipfun1='slice';
    end
    
    
    switch(lower(mipfun1))
        case 'mean'
            Xtmp=mean(Xtmp,2);
        case 'meanabs'
            Xtmp=mean(abs(Xtmp),2);
        case 'first'
            Xtmp=flip(Xtmp,2);
            [~,zi]=max(Xtmp>0,[],2);
            Xtmp=Xtmp(sub2ind(size(Xtmp),(1:size(Xtmp,1))',zi));
        case 'last'
            [~,zi]=max(Xtmp>0,[],2);
            Xtmp=Xtmp(sub2ind(size(Xtmp),(1:size(Xtmp,1))',zi));
        case 'slice'
            Xtmp=Xtmp(:,mipslice);
        otherwise
            [~,midx]=max(abs(Xtmp),[],2);
            Xtmp=Xtmp(sub2ind(size(Xtmp),1:size(Xtmp,1),midx'));
    end
    Xtmp=reshape(Xtmp,sz(setdiff(1:3,d)));
    mipimg{di}=flipud(permute(Xtmp,[2 1]));
end

dsz=cellfun(@size,mipimg,'uniformoutput',false);
maxsz=max(cat(1,dsz{:}),[],1);

mipimg=cellfun(@(x)padimageto(x,maxsz,[],bgcolor),mipimg,'uniformoutput',false);
