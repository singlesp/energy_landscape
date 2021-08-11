function w = tukey_outlier(x,k)
if(nargin == 1)
    k = 2;
end

rmad = median(x - min(x))/.6745;
rp = (x-min(x))/(k*rmad);
w = (1-rp.^2).^2;
