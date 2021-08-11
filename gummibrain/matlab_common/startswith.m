function s = startswith(strings, startstr)
s = [];

if(isempty(strings) || isempty(startstr))
    return;
end

N = numel(startstr);

if(ischar(strings))
    s = strncmp(strings,startstr,N);
    return;
end

if(iscell(strings))
    s = false(size(strings));
    strlen = cellfun(@numel,strings);
    validx = strlen >= N;
    
    s(validx) = cellfun(@(x)(all(x(1:N)==startstr)),strings(validx));
    
end