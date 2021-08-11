function displaycalib = ComputeCalibration(R,G,B,v,name)

if(nargin < 4 || isempty(v))
    v = linspace(0,1,numel(R));
end

if(nargin < 5 || isempty(name))
    name = '';
end

L0 = mean([R(1) G(1) B(1)]);
aR = R(end)-L0;
aG = G(end)-L0;
aB = B(end)-L0;

% 
% bidx = 3;
% bR = log((R(bidx)-L0)/(R(end)-L0))/log(v(bidx))
% bG = log((G(bidx)-L0)/(G(end)-L0))/log(v(bidx))
% bB = log((B(bidx)-L0)/(B(end)-L0))/log(v(bidx))


bR = log((R-L0)/(R(end)-L0))./log(v);
bG = log((G-L0)/(G(end)-L0))./log(v);
bB = log((B-L0)/(B(end)-L0))./log(v);

bR = mean(bR(2:end-1));
bG = mean(bG(2:end-1));
bB = mean(bB(2:end-1));

displaycalib.timestamp = datestr(now,'YYYYmmdd-HHMM');
displaycalib.text = name;
displaycalib.L0 = L0;
displaycalib.aR = aR;
displaycalib.aG = aG;
displaycalib.aB = aB;
displaycalib.bR = bR;
displaycalib.bG = bG;
displaycalib.bB = bB;
displaycalib.RGB = [R(:) G(:) B(:)];
displaycalib.values = v;
