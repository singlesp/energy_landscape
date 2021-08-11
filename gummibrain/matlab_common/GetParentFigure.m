function hfig = GetParentFigure(hObj)
hfig = hObj;

while(~strcmpi(get(hfig,'type'),'figure'))
    hfig = get(hfig,'Parent');
end
