function n = normrows(v)

n = v./repmat(sqrt(sum(v.^2,2)),1,3);