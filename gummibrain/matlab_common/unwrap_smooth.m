function a = unwrap_smooth(a,smoothfilt)
%if each column of "a" contains epochs that SHOULD have a steady phase, but it jumps between
% two (mult of 2*pi), then find one of those phases (a mode of the
% ksdensity), center and mod around that, then smooth it to get rid of the
% 2*pi jumps, then unwrap

badcol = round(std(a)/pi) >= 1;
h2b = a(:,badcol);
if(~isempty(h2b))
    for i = 1:size(h2b,2)
        [kx,xi] = ksdensity(h2b(:,i));
        [~,idx] = max(kx);
        hmode = xi(idx);
        h2b(:,i) = mod(h2b(:,i)-hmode+pi,2*pi)-pi + hmode;

    end
    h2b = filtfilt(smoothfilt,1,h2b);
    a(:,badcol) = h2b;
end
a = unwrap(a);