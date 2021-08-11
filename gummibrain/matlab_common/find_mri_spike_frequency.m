function freq = find_mri_spike_frequency(data,fs,freqrange)

if(nargin < 3)
    freqrange = [12 20];
end

F = abs(fft(data));
F = mean(F,2);
%freqs = [0:size(F,1)-1]*fs/size(F,1);
maxlag = ceil(max(freqrange)*size(F,1)/fs);


[xc lags] = xcorr(F,maxlag);

lagfreqs = lags*fs/size(F,1);
lag_range = lagfreqs >= min(freqrange) & lagfreqs <= max(freqrange);
[~,idx] = max(xc(lag_range));

lagfreq_range = lagfreqs(lag_range);
freq = lagfreq_range(idx);
