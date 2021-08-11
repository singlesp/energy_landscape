function wiener_kernel = wiener_deconv(x,y,wlength,wpost)
if(nargin < 4)
    wpost = 0;
end

wiener_size = wlength;

wiener_kernel = zeros(wiener_size,size(y,2));

for i = wiener_size-wpost+1:size(y,1)-wpost
    wiener_kernel = wiener_kernel + x(i+wpost-wiener_size:i+wpost-1,:).*repmat(y(i,:),wiener_size,1);
end
wiener_kernel = wiener_kernel./(size(y,1)-wiener_size);
