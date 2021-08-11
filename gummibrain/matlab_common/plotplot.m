function h = plotplot(x,y,params1,params2)
if(nargin < 4)
    params2 = {};
end


zmax = [];
ax = gca;
ch=get(ax,'children');
[h zmax] = plotaux(x,y,ch,zmax,params1);

if(~isempty(params2))
    np = get(ax,'nextplot');
    if(~strcmpi(np,'add'))
        set(ax,'nextplot','add');
    end
    [h2 zmax] = plotaux(x,y,ch,zmax,params2);
   
    h = [h; h2];
    if(~strcmpi(np,'add'))
        set(ax,'nextplot',np);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h zmax] = plotaux(x,y,ch,zmax,params)
charidx = cellfun(@ischar,params);

zidx = find(regexpimatch(params(charidx),'^z'));
if(~isempty(zidx))
    zval = params{charidx(zidx)+1};
    if(numel(zval) == 1)
        params{charidx(zidx)} = 'zdata';
        params{charidx(zidx)+1} = zval*ones(size(x));
    end
    zmax = max(zval);
elseif(isempty(zmax))
    zdata = get(ch,'zdata');
    if(isempty(zdata))
        zmax = 0;
    else
        z = cell2mat(cellfun(@max,zdata,'uniformoutput',false));
        zmax = max(z);
    end
end

aidx = regexpimatch(params(charidx),'alpha');
do_alpha = any(aidx);

if(do_alpha)
    h = plotalpha(x,y,params{:});
else
    h = plot(x,y,params{:});
end
if(~isempty(zmax))
    set(h,'zdata',2+zmax+0*get(h,'xdata'));
end