function V=nan2val(V,val)
if(~exist('val','var') || isempty(val))
    val=0;
end

V(~isfinite(V))=val;
