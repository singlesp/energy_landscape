% function snr = SpecSNR(F,snrsize,varargin)
%
% Compute spectral SNR (amp @ each freq divided by freqband around it)
%
% F = fft for all electrodes (each column is an electrode)
% snrsize = how wide (in indices) is the noise frequency band?
%
% optional params ('paramname',value): 
%   srate
%   snrsizeHz (to use this, leave snrsize blank and include the srate)

function varargout = SpecSNR(F,snrsize,varargin)
p = inputParser;
p.addParamValue('srate',[]);
p.addParamValue('snrsizehz',[]);
p.addParamValue('median',false);

p.parse(varargin{:});
r = p.Results;
do_median = r.median;

if(isempty(snrsize) && ~isempty(r.srate) && ~isempty(r.snrsizehz))
    freqbinsize = r.srate/size(F,1);
    snrsize = fix(.5 * r.snrsizehz / freqbinsize);
end

if(do_median)
    %F_noise = medfilt1(abs(F),2*snrsize+1,min(size(F,1),10000));
    F_noise = fastmedfilt1(abs(F),2*snrsize+1,min(size(F,1),10000));
else
    noisefilt = ones(2*snrsize+1,1);
    noisefilt(snrsize+1) = 0;
    noisefilt = noisefilt/sum(noisefilt);

    F_noise = filter(noisefilt,1,abs(F));
    F_noise = circshift(F_noise,[-snrsize 0]);
end
F_noise(round(end/2):end,:) = 1;
F_noise(F_noise == 0) = 1;
snr = F./abs(F_noise);

if(nargout > 1)
    varargout = {snr, snrsize};
else
    varargout = {snr};
end