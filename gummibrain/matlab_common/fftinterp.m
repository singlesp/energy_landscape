function fnew = fftinterp(F,freq0,varargin)
% fnew = fftinterp(F,freq0,['srate',<sampling freq>], ['window',<windowtype>]);
%
% F = complex fft
% freq0 = fft peak frequency (can be a vector of frequencies)
%   If srate is provided, assume this is a frequency in Hz
%   If no srate, assume this is an index (ie: Hz = (k-1)*fs/N )
% srate = sampling frequency in Hz
% window = type of temporal window used when first computing fft
%   window types: rectangular (default), hanning, hamming, blackman, blackman-harris (3 term)
%   window type determines which weightings for interpolation
%
% fnew = interpolated peak frequency
%   If srate provided, returns value in Hz
%   If no srate, return adjusted index (non-integer) s.t.: Hz = (k-1)*fs/N
% 
%
% from Jacobsen, E., and P. Kootsookos. 
% “Fast, Accurate Frequency Estimators [DSP Tips Tricks].” 
% IEEE Signal Processing Magazine 24, no. 3 (2007): 123–125. doi:10.1109/MSP.2007.361611.

p = inputParser;
p.addParamValue('window','');
p.addParamValue('srate',[]);

p.parse(varargin{:});
r = p.Results;

windowtype = r.window;
fs = r.srate;

N = size(F,1);

switch lower(windowtype)
    case {'','rect','rectangular'} %eq3 best
        P = [];
        Q = [];
    case {'hamming','hamm'} %eq4 best
        P = 1.22;
        Q = 0.6;
    case {'hanning','hann'} %eq4 best
        P = 1.36;
        Q = 0.55;
    case {'blackman','b'} %eq4 best
        P = 1.75;
        Q = 0.55;
    case {'blackmann-harris','bh','bh3'} %3 term
        P = 1.72;
        Q = 0.56;
end


%no starting freq provided.  Find overall FFT peak (in nyquist range)
if(nargin < 2 || isempty(freq0) || any(freq0 <= 0))
    [~,k] = max(abs(F(1:round(end/2))));
else
    if(isempty(fs) || fs <= 0)
        k = freq0;
    else
        k = round(1+freq0(:)/(fs/N));
        krange = -1:1;
        kvals = repmat(k,1,numel(krange))+repmat(krange,numel(k),1);
        [~,idx] = max(abs(F(kvals)),[],2);
        k = kvals(sub2ind(size(kvals),[1:numel(k)].',idx));
    end
end


%d = fftinterp(F,k0)
%fnew = fftinterp(F,fs,freq0)
%F, fs, k0
switch lower(windowtype)
    case {'','rect','rectangular'}
        %eq 3
        d3 = -real((F(k+1)-F(k-1))./(2*F(k)-F(k-1)-F(k+1)));
        d = d3;
    otherwise
        %eq 4
        d4 = P*(abs(F(k+1))-abs(F(k-1)))./(abs(F(k)) + abs(F(k-1)) + abs(F(k+1)));
        
        %eq5
        %d5 = real(Q*(F(k-1)-F(k+1))./(2*F(k)+F(k-1)+F(k+1)));
        d = d4;
end

k = reshape(k,size(freq0));
d = reshape(d,size(freq0));
if(isempty(fs) || fs <= 0)
    fnew = k+d;
else
    fnew = (k+d-1)*(fs/N);
end
