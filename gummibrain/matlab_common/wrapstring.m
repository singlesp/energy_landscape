function ws = wrapstring(s,width)

if(nargin < 2)
    width=80;
end

ws = repmat(' ',width,ceil(numel(s)/width));
ws(1:numel(s)) = s;

ws = ws';