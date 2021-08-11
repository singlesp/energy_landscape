clear all;
close all;
clc;

fs = 1000;

A_stop1 = 12;		% Attenuation in the first stopband
F_stop1 = 0.5;		% Edge of the stopband
F_pass1 = 2;	% Edge of the passband
F_pass2 = 30;	% Closing edge of the passband
F_stop2 = 35;	% Edge of the second stopband
A_stop2 = 12;		% Attenuation in the second stopband
A_pass = 1;		% Amount of ripple allowed in the passband
BandPassSpecObj = ...
    fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
    A_stop2, fs);
BandPassFilt = design(BandPassSpecObj,'butter');


LowPassSpecObj = ...
    fdesign.lowpass('Fp,Fst,Ap,Ast', ...
    4, 8, 1, 24,fs);

LowPassFilt = design(LowPassSpecObj,'butter');

freqz(BandPassFilt);


%%

close all;
trialsize = 30000;
N = 100;
Fsig_est = zeros(trialsize,1);
Ffilt_est = zeros(trialsize,1);
Fcplx_est = zeros(trialsize,1);


[filt_mag w] = freqz(BandPassFilt,trialsize/2);
filt_mag = [filt_mag; filt_mag(end:-1:1)];
figure;
freqs = linspace(0,fs,trialsize);
plot(freqs,10*log10(abs(filt_mag)));
return;
for i = 1:N
    signal = randn(trialsize,1);
    signal = signal-mean(signal);

    %signal = signal.*hann(numel(signal));
    %figure;
    %subplot(2,1,1);
    %plot(signal)
    signal0 = real(ifft(abs(fft(signal)))); %set all phases to 0
    signal0(1) = 0;
    %subplot(2,1,2);
    %plot(signal0);
    %return;
    signal = signal0;
    
    Fsig = fft(signal)/numel(signal);
    Fsig_est = Fsig_est + abs(Fsig);
    
    sigfilt = filter(BandPassFilt,signal);
    %sigfilt = filter(BandPassFilt,sigfilt(end:-1:1));
    %sigfilt = sigfilt(end:-1:1);
    
    Ffilt = fft(sigfilt)/numel(sigfilt);
    Ffilt_est = Ffilt_est + abs(Ffilt);
    
    Fcplx_est = Fcplx_est + angle(Ffilt);
end

freqs = linspace(0,fs,numel(signal));

Ffilt_ratio = abs(Ffilt_est./Fsig_est);

Fphs_est = -unwrap(Fcplx_est);

offset_est = fs*Fphs_est./(2*pi*freqs');
figure;
subplot(2,1,1);
plot(freqs,(Ffilt_ratio));
set(gca,'xlim',[0 30])
subplot(2,1,2);
%plot(freqs,angle(Fcplx_est));
%plot(freqs(1:end),-unwrap(Fcplx_est(1:end))./freqs(1:end)');
plot(freqs,offset_est.*Ffilt_ratio);
set(gca,'xlim',[0 30]);

