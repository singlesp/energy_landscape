function G = gausswindow(N,width,widthtype)

if(nargin < 2)
    width = N/2;
    widthtype = 'std';
elseif(nargin < 3)
    widthtype = 'fwhm';
end


if(strcmpi(widthtype,'fwhm'))
    gstd = width/(2*sqrt(2*log(2)));
else
    gstd = width;
end

x = 0:N-1;
gmean = x(end)/2;
G = exp(-((x-gmean).^2)/(2*gstd.^2));


G = G - min(G);
G = G/sum(G);
