function mx = spm_maxval(vol)
if isstruct(vol),
    mx = -Inf;
    for i=1:vol.dim(3),
        tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
        imx = max(tmp(isfinite(tmp)));
        if ~isempty(imx),mx = max(mx,imx);end
    end;
else
    mx = max(vol(isfinite(vol)));
end;