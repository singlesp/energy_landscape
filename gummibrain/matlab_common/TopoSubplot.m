function [handles valididx] = TopoSubplot(xvals, yvals, options)

handles = [];

if(~isfield(options,'locs'))
    return;
end

if(ischar(options.locs))
    L = readlocs(options.locs);
    options.locs = L;
end

if(nargin < 3)
    return;
end

if(~isfield(options,'axprop'))
    options.axprop = {};
end

if(~isfield(options,'plotprop'))
    options.plotprop = {};
end

if(~isfield(options,'masktype'))
    options.masktype = true;
end

if(~isfield(options,'showtitle'))
    options.showtitle = false;
end

if(~isfield(options,'showscale'))
    options.showscale = true;
end

fig = gcf;
locmask = StringMask({options.locs.labels},options.locmask,options.masktype);

handles = CreateTopoSubplots(options.locs,locmask,options.masktype);

for i = 1:numel(handles)
    if(~handles(i))
        continue;
    end

    set(fig,'CurrentAxes',handles(i));
    
    if(iscell(yvals))
        y = yvals{i};
    else
        y = yvals(:,i);
    end
    
    if(isempty(xvals))
        x = 1:numel(y);
    elseif(min(size(xvals)) == 1)
        x = xvals;
    else
        x = xvals(:,i);
    end
    plot(x,y,'hittest','off',options.plotprop{:});
    set(handles(i),'ButtonDownFcn', @clicker,options.axprop{:});
    
    %axis off;
    if(options.showtitle)
        title(options.locs(i).labels);
    end
end

if(nargout > 1)
    valididx = find(locmask);
end

if(options.showscale)
    axscale = copyobj(handles(1),gcf);
    pos = get(axscale,'position');
    set(axscale,'position',[.05 1-pos(4)-.05 pos(3:4)],'ButtonDownFcn',[]);
    hold(axscale,'on');
    cla(axscale);
    grid(axscale,'on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clicker(source, event)
newfig = figure;
plot(0,0);
pos = get(gca,'position');
delete(gca);

h = copyobj(source,newfig);
set(h,'position',pos,'ButtonDownFcn',[]);

ch = get(h,'children');
set(ch,'hittest','on');

