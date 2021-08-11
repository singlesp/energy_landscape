function s = uistruct(fig)

s = struct();

u = findobj(fig);
for i = 1:numel(u)
    tag = get(u(i),'tag');
    try
        s.(tag) = u(i);
    catch err
        % might have a tag name that can't be  used as a variable name
    end
end
