function [FPR,TPR,AUC,P]=fastroc(X,Y,dim)
class1=max(Y(:));
if(size(X,dim)~=numel(Y))
    error('Size of Y must match size(X,dim)');
end
[Yval,sortidx]=sort(X,dim);
Ysort=Y(sortidx)==class1;

TPR=1-bsxfun(@rdivide,cumsum(Ysort,dim),sum(Ysort,dim));
FPR=1-bsxfun(@rdivide,cumsum(1-Ysort,dim),sum(1-Ysort,dim));

dTPR=diff(TPR,[],dim);
if(dim==1)
    dTPR=dTPR([1 1:end],:);
    dTPR(1,:)=0;
elseif(dim==2)
    dTPR=dTPR(:,[1 1:end]);
    dTPR(:,1)=0;
else
    error('dim must be 1 or 2');
end
AUC=sum((FPR-1).*dTPR,dim);
P=nan(size(AUC));

%%
n1=(sum(Y==class1));
n2=(sum(Y~=class1));

%n1=min(n1,n2);
%n2=n1;

U=AUC*n1*n2;

umu=n1*n2/2;
usig=sqrt(n1*n2*(n1+n2+1)/12);

%this is not exact!
P=1-2*abs(normcdf(U,umu,usig)-.5);