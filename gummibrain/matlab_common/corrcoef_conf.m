function [rlo, rup] = corrcoef_conf(r, n, alpha)
% [rlo, rup] = corrcoef_conf(r, n, alpha)
% 
% Confidence interval for correlation coefficient
%
% r       = corrcoef (can be [m x 1] vector)
% n       = # of data points used to compute corrcoef
% alpha   = confidence interval to compute (default: .05) (can be [p x 1] vector)
% 
% rlo, rup = [m x p] matrix where each row contains CI for one R value

% (extracted from corrcoef.m)

if(nargin < 3 || isempty(alpha))
    alpha = .05;
end

rv = repmat(r(:),1,numel(alpha));
av = repmat(alpha(:),1,numel(r)).';

% Confidence bounds are degenerate if abs(r) = 1, NaN if r = NaN.
z = 0.5 * log((1+rv)./(1-rv));
zalpha = NaN(size(rv),class(r));
if(n <= 3)
    rlo = zalpha;
    rup = zalpha;
    return;
end

zalpha = (-erfinv(av - 1)) .* sqrt(2) ./ sqrt(n-3);

rlo = tanh(z-zalpha);
rup = tanh(z+zalpha);

