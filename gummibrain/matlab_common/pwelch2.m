function [pxx freqs] = pwelch2(x,varargin)

pxx = [];
freqs = [];
for i = 1:size(x,2)
    [p f] = pwelch(x(:,i),varargin{:});
    if(isempty(pxx))
        pxx = zeros(numel(p),size(x,2));
        freqs = f;
    end
    pxx(:,i) = p;
end