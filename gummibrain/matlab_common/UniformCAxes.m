function UniformCAxes(handles,color_symmetric,hsource)
if(nargin < 1)
    handles = gcf;
end

if(nargin < 2 || isempty(color_symmetric))
    color_symmetric = true;
end

if(nargin < 3)
    hsource = [];
end

axhandles = findobj(handles,'type','axes');
axtag = get(axhandles,'tag');
axleg = strcmpi(axtag,'Legend');
axcb = strcmpi(axtag,'Colorbar');
axsuptitle = strcmpi(axtag,'suptitle');
axhandles = axhandles(~axleg & ~axcb & ~axsuptitle);

if(isempty(hsource))
    hsource = axhandles;
end

%if(~isempty(findobj(handles,'type','image')))
if(~isempty(findobj(axhandles,'-property','CData')))
    %cl = cell2mat(get(hsource,'clim'));
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