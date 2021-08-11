function tilefigures(figs,mosaic)
if(~exist('figs','var') || isempty(figs))
    figs=findall(0,'type','figure');
    figs=sort(cell2mat(get(figs,'number')));
end
if(~exist('mosaic','var') || isempty(mosaic))
    n=numel(figs);
    n1=floor(sqrt(n));
    n2=ceil(n/n1);
    mosaic=[n1 n2];
end

screensz=get(0,'ScreenSize');
screensz=screensz(3:4);

figsz=[screensz(1)/mosaic(2) screensz(2)/mosaic(1)];

for i = 1:numel(figs)
    r=mosaic(1)-floor((i-1)/mosaic(2))-1;
    c=mod(i-1,mosaic(2));
    p=[figsz(1)*c figsz(2)*r figsz];
    set(figs(i),'outerposition',p);
    figure(figs(i));
end