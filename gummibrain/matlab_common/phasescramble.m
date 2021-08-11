function y = phasescramble(x, N)

if(nargin < 2)
    N = 1;
end

x = repmat(x,1,N);
F = fft(x);
Famp = abs(F);
Fphs = angle(F);
Fphs = reshape(shuffle(Fphs(:)),size(F));
y = real(ifft(Famp.*exp(sqrt(-1)*Fphs)));
y = y - repmat(mean(y,1),size(y,1),1) + repmat(mean(x,1),size(y,1),1);