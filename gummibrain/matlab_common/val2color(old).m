function colors = val2color(vals,arg2,c_lim)
if(nargin < 2)
    arg2 = gca;
end

if(numel(arg2) > 1) %it's a colormap array
    cmap = arg2;
    clim = [min(flatten(vals)) max(flatten(vals))];
    
elseif(strcmpi(get(arg2,'type'),'axes'))  %it's an axes object
    cmap = colormap;
    clim = get(arg2,'clim');
end

if(nargin > 2)
    clim = c_lim;
end

cmap_idx = min(size(cmap,1),max(1,fix((size(cmap,1)-1)*((vals - clim(1)) / (clim(2)-clim(1))) + 1)));
colors = reshape(cmap(cmap_idx,:),[size(vals) size(cmap,2)]);
