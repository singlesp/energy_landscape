function img = drawcolorbar(cmin,cmax,ticks,cmap,sz,varargin)
%img = drawcolorbar(cmin,cmax,ticks,cmap,sz,'param',value,...)
%
%Returns a rasterized colorbar image suitable for saving to file
%
%Inputs:
%   cmin
%   cmax
%   ticks
%   cmap
%   sz
%
%Options: 'paramname','value',...
%   tickformat
%   tickcolor
%   bgcolor
%   showbounds
%   fontsize

p = inputParser;
p.addParamValue('tickformat',[]);
p.addParamValue('tickcolor',[0 0 0]);
p.addParamValue('bgcolor',[1 1 1]);
p.addParamValue('showbounds',false);
p.addParamValue('fontsize',[]);
p.parse(varargin{:});
r = p.Results;
dumpstruct(r);


if(isempty(ticks))
    fig = figure('IntegerHandle','off','visible','off');
    hc=colorbar;
    set(gca,'clim',[cmin cmax]);
    set(hc,'units','pixels');
    p = get(hc,'position');
    p(4) = sz(1);
    set(hc,'position',p);
    ticks = get(hc,'ytick');
    close(fig);
end

if(showbounds)
    ticks = [ticks(:); cmin; cmax];
end
ticks = sort(unique(ticks));
    
fig = figure('IntegerHandle','off','color',bgcolor,'visible','off','position',[0 0 max(sz) 1.25*max(sz)]);
cmap_sz = size(cmap,1);
cimg = interp1(1:cmap_sz,cmap,linspace(1,cmap_sz,sz(1)));
cimg = reshape(cimg,[sz(1) 1 3]);
cimg = repmat(cimg, [1 sz(2) 1]);

tickpos = sz(1)*(ticks-cmin)/(cmax-cmin)+.5;
tickpos = min(max(tickpos,.5),sz(1)-.5);
if(~isempty(tickformat))
    ticks = cellfun(@(a)sprintf(tickformat,a),num2cell(ticks),'uniformoutput',false);
end

% toolboxerr = false;
% try
%     imshow(cimg);
% catch err
%     toolboxerr = true;
% end
toolboxerr=true;
if(toolboxerr)
    image(cimg);
    axis image;
end

ax = gca;

axis on;
axis xy;
set(ax,'ytick',tickpos,'yticklabel',ticks,'xtick',[],...
    'yaxislocation','right','tickdir','in','ycolor',tickcolor,'xcolor',tickcolor);
if(~isempty(fontsize))
    set(ax,'fontsize',fontsize);
end

if(true || toolboxerr)

    set(fig,'units','normalized');

    p=get(ax,'position');
    set(ax,'position',[p(1) .1 p(3) .78]);
end

img = export_fig(ax,'-a1');

%img = padimageto(CropBGColor(img,255*bgcolor),sz,5,255*bgcolor);
img = padimageto(img,sz,5,255*bgcolor);
close(fig);
