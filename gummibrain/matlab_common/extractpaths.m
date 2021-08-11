function paths = extractpaths(edges)
% paths = extractpaths(edges, open_ends)
%
% returns a connected sequence of vertex indices for a given set of
%   edges.  If there are multiple disconnected paths in this set of
%   edges, they will be listed sequentially, divided by a NaN.
%
% edges: Nx2 list of connected vertex indices

%find any points with only 1 connection (starting points for a non-loop path)
edges_u = sort(unique(edges(:)));
hc = histc(edges(:),edges_u);
open_ends = edges_u(hc==1);

if(isempty(open_ends))
    open_ends = edges(1);
end

cur_pos = open_ends(1);
open_ends = open_ends(2:end);

paths = cur_pos;

while(~isempty(edges))
    edgeidx = [find(edges(:,1) == cur_pos); find(edges(:,2) == cur_pos)];
    
    if(isempty(edgeidx))
        cur_is_end = open_ends == cur_pos;
        if(~isempty(cur_is_end))
            open_ends(cur_is_end) = [];
        end
        
        if(isempty(open_ends))
            if(isempty(edges))
                paths = [paths NaN];
                break;
            else
                cur_pos = edges(1);
                paths = [paths NaN cur_pos];
            end
        else
            cur_pos = open_ends(1);
            open_ends = open_ends(2:end);
            paths = [paths NaN cur_pos];
        end    
    else
        edgeidx = edgeidx(1); %only care about the first match
        cur_edge = edges(edgeidx,:);
        
        cur_pos = cur_edge(cur_edge ~= cur_pos);
        paths = [paths cur_pos];
        
        edges(edgeidx,:) = []; %remove this edge
    end
end

