function mn = spm_minval(vol)
if isstruct(vol),
    mn = Inf;
    for i=1:vol.dim(3),
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imn = min(tmp(isfinite(tmp)));
        if ~isempty(imn),mn = min(mn,imn);end
    end;
else
    mn = min(vol(isfinite(vol)));
end;