function n = normdim(v,d)

nvals = sqrt(sum(v.^2,d));
dimrep = ones(1,ndims(v));
dimrep(d) = size(v,d);

n = v./repmat(nvals,dimrep);
