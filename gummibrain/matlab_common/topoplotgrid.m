function h = topoplotgrid(topomat, varargin)
%topomat = ExT list where each column is a topo

cparam = StringIndex(varargin,'color');
if(~isempty(cparam) && cparam > 0)
    bgcolor = varargin{cparam+1};
    idx = true(size(varargin));
    idx([cparam cparam+1]) = false;
    varargin = {varargin{idx}};

    if(numel(bgcolor) == 1)
        bgcolor = bgcolor*[1 1 1];
    end
else
    bgcolor = .6*[1 1 1];
end

numtopo = size(topomat,2);
topoimg_grid = topoplotgridimage(topomat,varargin{:});

cols = ceil(sqrt(numtopo));
rows = ceil(numtopo/cols);
toposz = size(topoimg_grid,1)/rows;

yplus = strcat(num2str([0:cols:numtopo]'),'+');

im = imagesc(topoimg_grid);
set(im,'alphadata',1-isnan(topoimg_grid));
set(gca,'color',bgcolor);
set(gca,'xaxislocation','top','tickdir','out');
set(gca,'xtick',toposz*(1:cols)-toposz/2,'ytick',toposz*(1:rows)-toposz/2);
%set(gca,'xticklabel',1:cols,'yticklabel',1:rows);
%set(gca,'xticklabel',1:cols,'yticklabel',0:cols:numtopo);
set(gca,'xticklabel',1:cols,'yticklabel',yplus);
axis image;
UniformAxes(gcf,true);

h = gca;