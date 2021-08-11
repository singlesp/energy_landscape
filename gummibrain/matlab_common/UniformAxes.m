function UniformAxes(handles,color_symmetric,hsource,uniform_clim)

if(nargin < 1)
    handles = gcf;
end

if(~exist('color_symmetric','var') || isempty(color_symmetric))
    color_symmetric = false;
end

if(~exist('uniform_clim','var') || isempty(uniform_clim))
    uniform_clim = true;
end

if(~exist('hsource','var') || isempty(hsource))
    hsource = [];
end

axhandles = findobj(handles,'type','axes');

if(isempty(axhandles))
    return;
end
axtag = get(axhandles,'tag');
axleg = strcmpi(axtag,'Legend');
axcb = strcmpi(axtag,'Colorbar');
axsuptitle = strcmpi(axtag,'suptitle');
axhandles = axhandles(~axleg & ~axcb & ~axsuptitle);

if(isempty(hsource))
    hsource = axhandles;
end

xl = get(hsource,'xlim');
yl = get(hsource,'ylim');


if(numel(hsource) > 1)
    xl = cell2mat(xl);
    yl = cell2mat(yl);
end

set(axhandles,'xlim',[min(xl(:,1)) max(xl(:,2))]);
set(axhandles,'ylim',[min(yl(:,1)) max(yl(:,2))]);


if(uniform_clim && ~isempty(findobj(axhandles,'-property','CData')))
    cl = get(hsource,'clim');
    if(numel(hsource) > 1)
        cl = cell2mat(cl);
    end
    if(color_symmetric)
        cl = max(abs(cl(:)));
        cl = [-cl cl];
    else
        cl = [min(cl(:,1)) max(cl(:,2))];
    end
    set(axhandles,'clim',cl);
end
