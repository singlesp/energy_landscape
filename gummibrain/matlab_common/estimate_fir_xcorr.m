function h = estimate_fir_xcorr(x,y,N)
% sidex.m - Demonstration of the use of FFT cross-
% correlation to compute the impulse response 
% of a filter given its input and output.
% This is called "FIR system identification".
%
% adapted from : https://ccrma.stanford.edu/~jos/st/FIR_System_Identification.html

Nx = size(x,1); % input signal length 
%Nh = N; % filter length 
%Ny = Nx+Nh-1; % max output signal length
% FFT size to accommodate cross-correlation:
Nfft = 2^nextpow2(Nx); % FFT wants power of 2

%x = rand(1,Nx); % input signal = noise
%h = [1:Nh]; 	% the filter
xzp = [x; zeros(Nfft-size(x,1),size(x,2))]; % zero-padded input
%yzp = filter(h,1,xzp); % apply the filter
yzp = [y; zeros(Nfft-size(y,1),size(y,2))];

X = fft(xzp);   % input spectrum
Y = fft(yzp);   % output spectrum
Rxx = conj(X) .* X; % energy spectrum of x
Rxy = conj(X) .* Y; % cross-energy spectrum
Hxy = Rxy ./ Rxx;   % should be the freq. response
hxy = ifft(Hxy);    % should be the imp. response

h = hxy(1:N); 	    % print estimated impulse response
% freqz(hxy,1,Nfft);  % plot estimated freq response
% 
% err = norm(hxy - [h,zeros(1,Nfft-Nh)])/norm(h);
% disp(sprintf(['Impulse Response Error = ',...
% 	'%0.14f%%'],100*err));
% 
% err = norm(Hxy-fft([h,zeros(1,Nfft-Nh)]))/norm(h);
% disp(sprintf(['Frequency Response Error = ',...
% 	'%0.14f%%'],100*err));