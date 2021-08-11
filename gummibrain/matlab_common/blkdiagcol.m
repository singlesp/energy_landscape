function b = blkdiagcol(c)

b = cell(1,size(c,2));
for i = 1:size(c,2)
    b{i} = c(:,i);
end
b = blkdiag(b{:});