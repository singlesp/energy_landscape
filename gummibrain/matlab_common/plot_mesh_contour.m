function plot_mesh_contour(tri,verts,vals,levels,colors,linewidth)

if(nargin < 4 || isempty(levels))
    levels = median(vals);
end

if(nargin < 5 || isempty(colors))
    colors = [1 0 0];
end

if(nargin < 6)
    linewidth = 2;
end

hs = get(gca,'nextplot');
set(gca,'nextplot','add');

for l = 1:numel(levels)
    level = levels(l);
    color = colors(rem(l-1,size(colors,1))+1,:);
    
    [cx cy cz] = mesh_contour(tri,verts,vals,level);
    path_breaks = DivideNanVector(cx);
    num_paths = size(path_breaks,1);
    for i = 1:num_paths
        p = path_breaks(i,1):path_breaks(i,2);
        plot3(cx(p),cy(p),cz(p),'color',color,'linewidth',linewidth,'tag','contourline');
		%patch(cx(p),cy(p),cz(p),0,'edgecolor',color,'linewidth',linewidth,'edgealpha',.5,'facealpha',0,'tag','contourline');
    end

    %plot3(cx',cy',cz','color',color,'linewidth',linewidth,'tag','contourline');
end

set(gca,'nextplot',hs);