function varargout = plotstd(t,x,color,linestyle,fillstyle,bottomplot)

if(nargin == 1)
    x = t;
    t = 1:size(x,1);
end

t = t(:);
x_mean = mean(x,2);
%x_ste = ste(x,2);
x_ste = std(x,[],2);

if(nargin < 3)
    color = 'b';
end

if(nargin < 4 || isempty(linestyle))
    linestyle = {};
end

if(nargin < 5 || isempty(fillstyle))
    fillstyle = {};
end

if(nargin < 6 || isempty(bottomplot))
    bottomplot = false;
end

linestyle = ['linewidth', 2, linestyle];
fillstyle = ['edgealpha', 0,'facealpha', .5, fillstyle];


do_hold = get(gca,'nextplot');

if(bottomplot)
    %x_min = min(x_mean);
    x_min = 0;
    hfill = fill([t; t([end 1])],[x_min+x_ste; x_min; x_min],color,fillstyle{:});
else
    hfill = fill([t; t(end:-1:1)],[x_mean+x_ste; x_mean(end:-1:1)-x_ste(end:-1:1)],color,fillstyle{:});
end
hold on;
hplot = plot(t,x_mean,'color',color,linestyle{:});

set(gca,'nextplot',do_hold);

if(nargout == 1)
	varargout = {hplot};
elseif(nargout == 2)
	varargout = {hplot,hfill};
end