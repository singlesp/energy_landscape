function [R P RLO RUP lags] = xcorr_stats(varargin)
a =StringIndex(varargin,'alpha');
if(~isempty(a))
    alpha = varargin{a+1};
    varargin(a:a+1) = [];
else
    alpha = .05;
end
x = varargin{1};
n = size(x,1);

[r l] = xcorr(varargin{:});
Tstat = r .* sqrt((n-2) ./ (1 - r.^2));
p = 2*tpvalue(-abs(Tstat),n-2);
p(l==0) = r(l==0);

z = 0.5 * log((1+r)./(1-r));
z(l==0) = 0;
zalpha = NaN(size(n),class(x));
if any(n>3)
    zalpha(n>3) = (-erfinv(alpha - 1)) .* sqrt(2) ./ sqrt(n(n>3)-3);
end
rlo = tanh(z-zalpha);
rup = tanh(z+zalpha);
rlo(l==0) = r(l==0);
rup(l==0) = r(l==0);

R = r;

if(nargout > 1)
    P = p;
end
if(nargout > 3)
    RLO = rlo;
    RUP = rup;
end

if(nargout > 4)
    lags = l;
end


%%
function [p rlo rup] = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic.

% copied from matlab corrcoef();

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;

if(nargout > 1)
    
end