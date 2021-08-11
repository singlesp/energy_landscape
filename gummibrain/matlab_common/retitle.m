function varargout = retitle(figs, show_number)

if(nargin < 1 || isempty(figs))
    figs = gcf;
end

if(strcmpi(figs,'all'))
    figs = findobj(0,'type','figure');
end

if(nargin < 2)
    show_number = false;
end

titlestr = [];
for f = 1:numel(figs)
    titlestr{f} = '';
    fig = figs(f);
    
    st=findobj(fig,'tag','suptitle');
    if(~isempty(st))
        titlestr{f} = get(findobj(st,'type','text'),'string');
    else
        children = findobj(fig,'type','axes');
        titlestr{f} = get(get(children(end),'title'),'string');
    end
    if(~isempty(titlestr{f}) && iscell(titlestr{f}))
        titlestr{f} = titlestr{f}{1};
    end

end

numtitle = {'off','on'};
for f = 1:numel(figs)
    if(isempty(titlestr{f}))
        continue;
    end
    
    set(figs(f),'NumberTitle',numtitle{show_number+1});
    set(figs(f),'Name',titlestr{f});
end

if(nargout > 0)
    varargout = figs;
end
