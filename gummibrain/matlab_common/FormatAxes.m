function FormatAxes(handles, varargin)

axhandles = findobj(handles,'type','axes');
if(~ischar(varargin{1}) || strcmpi(varargin{1},'off'))
    axis(axhandles,'off');
    varargin = {varargin{2:end}};
else
    set(axhandles,varargin{:});
end

