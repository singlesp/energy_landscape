function mask = StringMask(allnames,names,masktype)


if(nargin < 3 || isempty(masktype))
    masktype = true;
end

mask = masktype*ones(size(allnames));

if(nargin < 2 || isempty(names))
    mask(:) = true;
    return;
end

%%% if locsmask was an array of electrode names, DISPLAY those whos names
%       are given.  if first item is 'bad', REMOVE the given electrodes 
if(iscell(names))

    idx = StringIndex(allnames,names);
    
    mask = false(size(allnames));
    mask(idx) = true;
    
    if(~masktype)
        mask = ~mask;
    end
end
