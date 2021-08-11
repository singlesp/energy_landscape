function F = idft(y, freq, N, srate, dim)
%
%Evaluate CONTINUOUS (freq) Fourier transform at specific frequency
%
%F = idft(y, freq, [srate])
%
%y = TxChannels signal
%freq = 1xN frequencies to evaluate
%srate = sampling freq.  If not specified, assume freq is already 
%   in sample radians (2*pi*freq/srate)
%
%returns F = NxChannels 

if(nargin >= 4 && ~isempty(srate) && srate > 0)
    w = 2*pi*freq/srate;
else
    w = freq;
end

if(nargin < 5 || isempty(dim))
    dim = 1;
end

%%

sz = size(y);
otherdims = setdiff(1:numel(sz),dim);
dimperm1 = [dim otherdims];
dimperm2 = zeros(size(dimperm1));
dimperm2([dim otherdims]) = 1:numel(dimperm1);

y = reshape(permute(y,dimperm1),sz(dim),[]);

%%
%% only diff from cft is negative exponent and divide by N
n=[0:N-1];
%F=exp(1i*w(:)*n)*y/size(y,1);
F = exp(1i*w(:)*n)*y/N + exp(1i*(2*pi-w(:))*n)*conj(y)/N;

%F2 = 2*(cos(w(:)*n)*real(y)/N - sin(w(:)*n)*imag(y)/N);

%%

F = permute(reshape(F,[N sz(otherdims)]),dimperm2);
