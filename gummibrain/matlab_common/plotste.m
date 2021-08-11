function varargout = plotste(t,x,color,linestyle,fillstyle)

if(nargin == 1)
    x = t;
    t = 1:size(x,1);
end

t = t(:);
x_mean = mean(x,2);
x_ste = ste(x,2);
%x_ste = std(x,[],2);

if(nargin < 3)
    color = 'b';
end

if(nargin < 4)
    linestyle = {};
end

if(nargin < 5)
    fillstyle = {};
end

linestyle = ['linewidth', 2, linestyle];
fillstyle = ['edgealpha', 0,'facealpha', .5, fillstyle];


do_hold = get(gca,'nextplot');

hfill = fill([t; t(end:-1:1)],[x_mean+x_ste; x_mean(end:-1:1)-x_ste(end:-1:1)],color,fillstyle{:});
hold on;
hplot = plot(t,x_mean,'color',color,linestyle{:});

set(gca,'nextplot',do_hold);

if(nargout == 1)
	varargout = {hplot};
elseif(nargout == 2)
	varargout = {hplot,hfill};
end