function [eigvec, sv] = minipca(data)
%function [eigvec sv] = minipca(data)
%
% KJ adapted from MathWorks pca.m

%center and decompose
data=bsxfun(@minus,data,nanmean(data,1));
[eigvec, eigValueDiag] = eig(cov(data));

%sort by eigen value
[sv, idx] = sort(diag(eigValueDiag), 'descend');
eigvec = eigvec(:, idx);

% Enforce a sign convention on the coefficients -- the largest element in
% each column will have a positive sign.
[~,maxind] = max(abs(eigvec), [], 1);
[d1, d2] = size(eigvec);
colsign = sign(eigvec(maxind + (0:d1:(d2-1)*d1)));
eigvec = bsxfun(@times, eigvec, colsign);
