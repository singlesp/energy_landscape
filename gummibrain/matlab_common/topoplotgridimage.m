function img = topoplotgridimage(topomat, varargin)
%topomat = ExT list where each column is a topo


topoimg_all = [];
numtopo = size(topomat,2);

padsz = 1;

for t = 1:numtopo
    [h timg] = topoplot(topomat(:,t),varargin{:},'noplot','on');
    topoimg_all(:,:,t) = padarray(timg,[padsz padsz],nan);
end

toposz = size(topoimg_all,1);
cols = ceil(sqrt(numtopo));
rows = ceil(numtopo/cols);
topoimg_grid = nan(toposz*rows,toposz*cols);

for t = 1:numtopo
    r = floor((t-1)/cols);
    c = mod(t-1,cols);
    topoimg_grid((toposz*r:toposz*(r+1)-1)+1,(toposz*c:toposz*(c+1)-1)+1) = topoimg_all(end:-1:1,:,t);
end

img = padarray(topoimg_grid,[padsz padsz],nan);