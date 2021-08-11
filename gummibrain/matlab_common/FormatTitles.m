function FormatTitles(handles,varargin)

ax = findobj(handles,'type','axes');
args = varargin;

if(~ischar(varargin{1}) || strcmpi(varargin{1},'off'))
    args = {'string',[]};
end

cellfun(@(t)(set(t,args{:})),get(ax,'title'));
