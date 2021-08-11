function clearchildren(hObj)
ch = get(hObj,'children');
if(isempty(ch))
    return;
end
v = get(ch,'handlevisibility');
if(ischar(v))
    v = {v};
end
i = StringIndex(v,'on');
delete(ch(i));
