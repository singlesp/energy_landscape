function X = submean(X,dim)
if(nargin < 2 || isempty(dim))
    dim = 1;
end
X = bsxfun(@minus,X,mean(X,dim));