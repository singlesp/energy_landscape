function [G gx] = spherespline(Xh,Yh,Zh,Xe,Ye,Ze)
% taken from eeglab headplot

fprintf('Setting up splining matrix.\n');
enum = length(Xe);
onemat = ones(enum,1);
G = zeros(enum,enum);
for i = 1:enum
    ei = onemat-sqrt((Xe(i)*onemat-Xe).^2 + (Ye(i)*onemat-Ye).^2 + ...
                     (Ze(i)*onemat-Ze).^2); % default was /2 and no sqrt
    gxt = zeros(1,enum);
    for j = 1:enum
        gxt(j) = calcgx(ei(j));
    end
    G(i,:) = gxt;
end
   
fprintf('Calculating splining matrix.\n');

gx = fastcalcgx(Xh,Yh,Zh,Xe,Ye,Ze);
        
%%%%%%%%%%%%%%%%%%%
function [out] = calcgx(in)

out = 0;
m = 4;       % 4th degree Legendre polynomial
for n = 1:7  % compute 7 terms
  L = legendre(n,in);
    out = out + ((2*n+1)/(n^m*(n+1)^m))*L(1);
end
out = out/(4*pi);

    
%%%%%%%%%%%%%%%%%%%
function gx = fastcalcgx(x,y,z,Xe,Ye,Ze)

onemat = ones(length(x),length(Xe));
EI = onemat - sqrt((repmat(x,1,length(Xe)) - repmat(Xe',length(x),1)).^2 +... 
                    (repmat(y,1,length(Xe)) - repmat(Ye',length(x),1)).^2 +...
                    (repmat(z,1,length(Xe)) - repmat(Ze',length(x),1)).^2);
%
EI = EI./max(abs(EI(:)));
gx = zeros(length(x),length(Xe));
m = 4;
%hwb = waitbar(0,'Computing spline file (only done once)...');
hwbend = 7;
for n = 1:7
    L = legendre(n,EI);
    gx = gx + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
    %waitbar(n/hwbend,hwb);
end
gx = gx/(4*pi);    
%close(hwb);