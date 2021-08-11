function S = catstruct(Sarray)

if(~isstruct(Sarray))
    error('inputs must be struct array');
end

if(numel(Sarray) <= 1)
    S = Sarray;
    return;
end

S = struct();
fn = fieldnames(Sarray(1));
for f = 1:numel(fn)
    newval = [];
    is_newval = false;
    sval = {Sarray.(fn{f})};
    %ctype = cellfun(@class,sval,'uniformoutput',false);

    
    if(all(cellfun(@isnumeric,sval)))
        sval = sval(~cellfun(@isempty,sval));
        if(isempty(sval))
            newval = [];
            is_newval = true;
        else
            nd = cellfun(@ndims,sval);
            if(all(nd == nd(1)))
                nd = nd(1);
                csz = cellfun(@size,sval,'uniformoutput',false);
                sz = cat(1,csz{:});
                eqdim = false(nd,1);
                for d = 1:nd
                    eqdim(d) = all(sz(:,d) == sz(1,d));
                end
                catdim = find(~eqdim);
                if(numel(catdim) == 0)
                    newval = cat(1,sval{:});
                    is_newval = true;
                elseif(numel(catdim == 1))
                    newval = cat(catdim,sval{:});
                    is_newval = true;
                end
            end
        end
    end
    
    if(~is_newval)
        newval = sval;
    end
    S.(fn{f}) = newval;
end