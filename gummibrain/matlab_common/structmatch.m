function [Smatch, idx] = structmatch(structlist, searchstruct, strcmpfcn, match_true_or_false)

if(isempty(structlist))
    Smatch = [];
    if(nargout > 1)
        idx = [];
    end
    return;
end

if(nargin < 3 || isempty(strcmpfcn))
    strcmpfcn = @strcmp;
end

if(ischar(strcmpfcn))
    strcmpfcn = eval(['@' strcmpfcn]);
end

if(nargin < 4 || isempty(match_true_or_false))
    match_true_or_false = true;
end

fields = fieldnames(searchstruct);
m = true(size(structlist));
for f = 1:numel(fields)
    fn = fields{f};
    v = searchstruct.(fn);
    if(isempty(v))
        r = cellfun(@isempty,{structlist.(fn)});
    elseif(isnumeric(v) || islogical(v))
        r = [structlist.(fn)] == v;
    elseif(ischar(v))
        r = strcmpfcn({structlist.(fn)},v);
    else
        error('cannot compare field of type: %s',class(v));
    end
    m = m & reshape(r,size(m));
end
if(~match_true_or_false)
    m=~m;
end

Smatch = structlist(m);
if(nargout > 1)
    idx = find(m);
end