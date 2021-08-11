function s_out = scalarstring(str,b,forcesign)

if(nargin < 3 || isempty(forcesign))
    forcesign = false;
end

s_out = {};
for i = 1:numel(b)
    if(b(i) == 0)
        s = '';
    elseif(b(i) == 1)
        s = str;
    elseif(b(i) == -1)
        s = strcat('-',str);
    else
        s = sprintf('%d%s',b(i),str);
    end

    if(forcesign && b(i) > 0)
        s = strcat('+',s);
    end
    s_out{i} = s;
end

if(numel(b) == 1)
    s_out = s_out{1};
end