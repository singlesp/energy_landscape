function C = fastmatrixcorr(data,do_demean)
if(~exist('do_demean','var') || isempty(do_demean))
    do_demean=false;
end
if(do_demean)
    data=demean(data,2);
end
dn=sqrt(sum(data.^2,2));
dn(dn<1e-6)=median(dn);
datanorm=bsxfun(@rdivide,data,dn);
C=datanorm*datanorm';
