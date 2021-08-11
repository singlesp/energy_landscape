function axclick_leftright(src,ev,callback_functions)
p = get(src,'CurrentPoint');
p = p(1,1:2);
numel(callback_functions)
if(numel(callback_functions) > 1 && strcmpi(get(gcbf,'selectiontype'),'alt')) %right click
    callback_functions{2}(p);
else
    callback_functions{1}(p);
end