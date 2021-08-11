function F = dft(y, freq, srate,dim)
%
%Evaluate CONTINUOUS (freq) Fourier transform at specific frequency
%
%F = dft(y, freq, [srate])
%
%y = TxChannels signal
%freq = 1xN frequencies to evaluate
%srate = sampling freq.  If not specified, assume freq is already 
%   in sample radians (2*pi*freq/srate)
%
%returns F = NxChannels 

if(nargin >= 3 && ~isempty(srate) && srate > 0)
    w = 2*pi*freq/srate;
else
    w = freq;
end

if(nargin < 4 || isempty(dim))
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

n=[0:size(y,1)-1];
F=exp(-1i*w(:)*n)*y;

%F1=exp(-1i*(w(:))*n)*y;
%F2=exp(-1i*(2*pi-w(:))*n)*conj(y);
%%

F = permute(reshape(F,[numel(w) sz(otherdims)]),dimperm2);
