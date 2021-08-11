function [TFR,timeVec,freqVec] = morletTFR(S,freqVec,Fs,width)
S = S';
 

timeVec = (1:size(S,2))/Fs;  

B = zeros(length(freqVec),size(S,2)); 

for i=1:size(S,1)          
     %fprintf(1,'%d ',i);
     %if(~mod(i,20))
     %    fprintf('\n');
     %end
    for j=1:length(freqVec)
       B(j,:) = energyvec(freqVec(j),detrend(S(i,:)),Fs,width) + B(j,:);
    end
end
TFR = B/size(S,1);     


function y = energyvec(f,s,Fs,width)
dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = conv(s,m);
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));



function y = morlet(f,t,width)
sf = f/width;
st = 1/(2*pi*sf);
A = 1/(st*sqrt(2*pi));

y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);
