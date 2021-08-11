function [ax freqs] = fftplot(F,fs,freqrange,varargin)

if(nargin < 3 || isempty(freqrange))
    freqrange = [0 30];
end

if(numel(freqrange) == 1)
    freqrange = [0 freqrange];
end

if(any(size(F)==1))
    F = F(:);
end

freqs = [0:size(F,1)-1]*fs/size(F,1);
stem(freqs,F,'marker','none',varargin{:});
ax = gca;
set(ax,'xlim',freqrange);