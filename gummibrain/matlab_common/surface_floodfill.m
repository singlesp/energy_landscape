%surf_colors = surface_floodfill(seedvert_idx, verts,tris,color1,color2)
% replace color1 with color2
function new_colors = surface_floodfill(seedvert_idx, verts, tris, vert_colors, color1, color2, max_dist)

    if(nargin < 7 || max_dist < 0)
        max_dist_sq = -1;
    else
        max_dist_sq = max_dist*max_dist;
    end

    color_verts = seedvert_idx;
    
    while(numel(color_verts) > 0)
        new_color_verts = [];
        for i = 1:numel(color_verts)
            cv = color_verts(i);
            
            vert_colors(cv) = color2;
            
            v = connected_verts(cv,tris); %find all verts connected to seed
            if(~isempty(color1))
                v = v(vert_colors(v) == color1); %only connected verts with right color (if fill mode)
            else
                v = v(vert_colors(v) ~= color2); %only connected verts not yet filled in (if paint mode)
            end
            
            if(max_dist_sq > 0) %only connected verts within distance from original seed
                d2 = sum((verts(v,:)-ones(numel(v),1)*verts(seedvert_idx,:)).^2,2);
                v = v(d2 < max_dist_sq);
            end

            new_color_verts = [new_color_verts; v];
        end
        color_verts = unique(new_color_verts);
    end
    
    new_colors = vert_colors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = connected_verts(seedvert_idx,tris)

    new_tris = tris(any(tris == seedvert_idx,2),:);
    v = unique(new_tris(:));