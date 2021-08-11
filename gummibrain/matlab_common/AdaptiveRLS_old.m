function [Hcn,Hsn] = AdaptiveRLS(pin,y,f0,order,fs)
%adaptive RLS filter

% lamda =1;
% f0 = 10;
% fs = 1000;
% order = 2;
% N = 1000;
% y = 0.5*sin(2*pi*order*f0*(0:N-1)/fs+pi/6)+2*rand(1,N)-1;
% pin = lamda/(1-lamda);
lamda = 1;%pin/(pin+1);
N = numel(y);
% N = N-mod(N,pin);
k = 0;
for n = pin:N-1
%for n = pin:N-1
    k = k+1;
    if n-pin<0
        startN = 0;
    else
        startN = n-pin+1;
    end

    w_idx = startN:n;
    w_lamda = lamda.^(n-w_idx);
    %w_lamda = hanning(numel(w_idx))';

    w_t = 2*pi*order*f0*w_idx/fs;
    w_sin = sin(w_t);
    w_cos = cos(w_t);

    %Hcn(k) = sum(w_lamda.*y(w_idx+1).*w_sin) / sum(w_lamda.*(w_sin.^2));
    %Hsn(k) = sum(w_lamda.*y(w_idx+1).*w_cos) / sum(w_lamda.*(w_cos.^2));

    %%%%%%%%%%%%%% why +1 for signal indices vs reference sin/cos?
    %%%%%%%%%%%%%% why using lambda = 1?
    %%%%%%%%%%%%%% why using epoch size = 1200, not shorter?
    %%%%%%%%%%%%%% 
    Hcn(k) = ...
       sum((lamda.^(n-(startN:n)).*y(startN+1:n+1)).*sin(2*pi*order*f0*(startN:n)/fs))./...
       sum((lamda.^(n-(startN:n)).*(sin(2*pi*order*f0*(startN:n)/fs))).^2);
    Hsn(k) = ...
       sum((lamda.^(n-(startN:n)).*y(startN+1:n+1)).*cos(2*pi*order*f0*(startN:n)/fs))./...
       sum((lamda.^(n-(startN:n)).*(cos(2*pi*order*f0*(startN:n)/fs))).^2);
end