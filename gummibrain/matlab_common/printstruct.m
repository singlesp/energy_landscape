function text = printstruct(S)

fnames = fieldnames(S);
text = {};
for f = 1:numel(fnames)
    fname = fnames{f};
    val = S.(fname);
    valtxt = printvariable(val);
    text{f} = [fname ' = ' valtxt ';'];
end

function txt = printnumber(v)
isflipped = false;
if(size(v,1) > 1 && size(v,2) > 1)
    warning('skipping... cant print 2D matrices yet');
    txt = sprintf('<%dx%d matrix>',size(v));
    return;
elseif(size(v,1) > 1)
    warning('flipping dimension for printing');
    isflipped = true;
    v = v';
end
txt = num2str(v);
%if(numel(v) > 1)
    txt = ['[' txt ']'];
    if(isflipped)
        txt = [txt ''''];
    end
%end

function txt = printstring(v)
txt = ['''' v ''''];

function txt = printlogical(v)
s = {'false','true'};
txt = s(v+1);
txt(1:end-1) = cellfun(@(x)([x ' ']), txt(1:end-1),'uniformoutput',false);
txt = [txt{:}];
if(numel(v) > 1)
    txt = ['[' txt ']'];
end

function txt = printcell(c)
txt = '{';
for f = 1:numel(c)
    txt = [txt printvariable(c{f})];
    if(f < numel(c))
        txt = [txt ', '];
    end
end
txt = [txt '}'];

function txt = printvariable(v,maxarray)

if(nargin < 2)
    maxarray = 30;
end

%txt = '<cannot display>';


if(~ischar(v) && (numel(v) > maxarray || (size(v,1) > 1 && size(v,2) > 1)))
    txt = sprintf('<%dx%d %s>',size(v), class(v));
else
    
    if(ischar(v))
        txt = printstring(v);
    elseif(isnumeric(v))
        txt = printnumber(v);
    elseif(islogical(v))
        txt = printlogical(v);
    elseif(iscell(v))
        txt = printcell(v);
    %elseif(isstruct(v))
    %    txt = printstruct(v);
    else
        txt = sprintf('<%dx%d %s>',size(v), class(v));
    end
end



