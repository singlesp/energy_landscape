function varargout = toptitle(h,varargin)

if(nargin < 1)
    h = gcf;
end

%if(numel(varargin) == 1 && iscell(varargin{1}) && numel(varargin{1}) == 1)
if(numel(varargin) == 1 && iscell(varargin{1}))
    varargin = varargin{1};
end

if(numel(varargin) == 1 && ischar(varargin{1}) && strcmpi(varargin{1},'inner'))
    params = {'verticalalignment','top'};
elseif(numel(varargin) == 1 && ischar(varargin{1}) && strcmpi(varargin{1},'middle'))
    params = {'verticalalignment','middle'};
elseif(numel(varargin) == 1 && ischar(varargin{1}) && strcmpi(varargin{1},'outer'))
    params = {'verticalalignment','bottom'};
else
    params = varargin;
end

if(iscell(h))
    h = cell2mat(h);
end

ht_all = [];
for i = 1:numel(h)
    if(~ishandle(h(i)))
        continue;
    end
    ax = [];
    ht = [];
    hfig = [];
    if(strcmpi(get(h(i),'type'),'figure'))
        ht = toptitle(gaa(h(i)),params{:});
    elseif(strcmpi(get(h(i),'type'),'axes'))
        ax = h(i);
        ht = get(ax,'title');
    elseif(strcmpi(get(h(i),'type'),'text'))
        ax = get(h(i),'parent');
        ht = h(i);
    end
    
    
    if(isempty(ax) || isempty(ht))
        continue;
    end
    
    ht_all = [ht_all ht];
    p = get(ht,'position');
    p(2) = max(get(ax,'ylim'));
    set(ht,'position',p,params{:});

    hfig = ax;
    while get(hfig,'parent') ~= 0
        hfig = get(hfig,'parent');
    end
    
    ttparam = getappdata(hfig,'toptitle');
    if(isempty(ttparam))
        setappdata(hfig,'toptitle',params);
        iptaddcallback(hfig,'ResizeFcn',@(src,ev)(toptitle(src,getappdata(src,'toptitle'))));
        %iptaddcallback(hfig,'ResizeFcn',@(src,ev)(toptitle(src,'inner')));
    end
end

if(nargout > 0)
    varargout = {ht_all};
end