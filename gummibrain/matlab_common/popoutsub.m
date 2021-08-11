function varargout = popoutsub(ax)
if(nargin < 1)
    ax = gca;
end

hfig = figure;
newax = copyobj(ax,hfig);
set(newax,'position',[0.1 0.1 .8 .8]);
%if(~isempty(get(ax,'title')))

if(nargout == 1)
    varargout = {newax};
end