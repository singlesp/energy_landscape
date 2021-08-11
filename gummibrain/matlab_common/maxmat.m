function varargout = maxmat(A)

[~,idx] = max(A(:));
[r1 r2 r3 r4 r5 r6] = ind2sub(size(A),idx);
r = [r1 r2 r3 r4 r5 r6];

if(nargout == 1)
    varargout = {r(1:ndims(A))};
else
    varargout = num2cell(r(1:nargout));
end