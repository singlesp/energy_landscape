function hsurf = clean_headplot(ax, show_electrodes, color_blend_mode, cminmax, reset, cmap)
%hsurf = clean_headplot(ax, show_electrodes, color_blend_mode, cminmax, reset, cmap)

if(nargin < 1)
    ax = gca;
end

%if preevious stored blending data for this headplot, don't have to
%recalculate everything again (time consuming)

has_userdata = ~isempty(get(ax,'userdata'));
if(has_userdata && ~reset)
    userdata = get(ax,'userdata');
else
    userdata = struct();
end

if(nargin < 2)
    show_electrodes = true;
end

if(nargin < 3 || isempty(color_blend_mode))
    color_blend_mode = 'blend';
end
color_blend_mode = lower(color_blend_mode);

if(nargin < 4)
    cminmax = [];
end

if(nargin < 6)
    cmap = [];
end

ch = get(ax,'children');
hsurf = ch(strcmpi(get(ch,'type'),'patch'));
hline = ch(strcmpi(get(ch,'type'),'line'));
hcontour = hline(strcmpi(get(hline,'tag'),'contourline')); %get contour objects
hline = hline(strcmpi(get(hline,'tag'),'')); %ignore contour objects

htext = ch(strcmpi(get(ch,'type'),'text'));


surf_loc = get(hsurf,'Vertices');
surf_tri = get(hsurf,'Faces');
surf_val = get(hsurf,'FaceVertexCData');

%%%%%%%%%%%%%%%%%
num_surfloc = size(surf_loc,1);

if(has_userdata)
    blend_alpha = userdata.blend_alpha;
else
    if(~isempty(htext))
        %eeg_loc = cell2mat(get(htext,'Position'));
        eeg_loc = cell2mat([get(hline(:),'XData'); get(hline(:),'YData'); get(hline(:),'ZData')]);
        eeg_loc = reshape(eeg_loc(:,1),size(eeg_loc,1)/3,3);
    else
        eeg_loc = [get(hline(1),'XData'); get(hline(1),'YData'); get(hline(1),'ZData')]';
    end

    %find scalp verts that fall within the convex hull defined by electrodes
    khull = convhulln(eeg_loc);
    surf_incap = inhull(surf_loc,eeg_loc,khull);

    %find scalp triangles that fall within this hull
    tri_incap = inhull(surf_loc(surf_tri(:),:),eeg_loc,khull);
    tri_incap = reshape(tri_incap,size(surf_tri));
    tri_incap = all(tri_incap,2);

    %find the open edges of this set of faces (the edge of the cap scalp surface)
    openedge = surfedge(surf_tri(tri_incap,:));
    loops = extractloops(openedge); %breaks into connected loops

    loop_breaks = DivideNanVector(loops);
    num_loops = size(loop_breaks,1);

    cap_edge = [];
    for i = 1:num_loops
        loop = loops(loop_breaks(i,1):loop_breaks(i,2));
        cap_edge = [cap_edge; surf_loc(loop,:)];
    end

    %blend away the scalp coloration as you get outside of the electrode cap
    % (based on the distance of scalp verts from the cap's edge defined above
    cap_dist = distance(surf_loc',cap_edge');
    sort_dist = sort(cap_dist,2);

    blend_width = 10;
    blend_center = blend_width;
    blend_alpha = 1./(1+exp((sort_dist(:,1)-blend_center)/blend_width));

    blend_alpha(surf_incap) = 1; %make sure everything INSIDE the cap is colored

    userdata.blend_alpha = blend_alpha;
    switch(color_blend_mode)
        case 'erase'
            blend_alpha(:) = 0;
        case {'noblend','none'}
            blend_alpha(:) = 1;
    end
end

%already have RGB values, just do blend
if(size(surf_val,2) == 3) 
    C = surf_val;
%just have data values, need to create a colormap manually so we can control
%blending, etc.    
else
    if(isempty(cmap))
        cmap_size = 256;
        %cmap = redgrayblue(cmap_size);
        cmap = jet(cmap_size);
    else
        cmap_size = size(cmap,1);
    end

    if(isempty(cminmax))
        m = max(abs(surf_val));
        surf_min = -m;
        surf_max = m;
    else
        surf_min = min(cminmax);
        surf_max = max(cminmax);
    end
    
    surf_range = surf_max-surf_min;

    cval = fix((cmap_size-1)*(surf_val - surf_min)/surf_range)+1;
    cval = max(min(cval,cmap_size),1);

    C = cmap(cval,:);
end

%(have to shape color matrix in a weird way)
C = reshape(C,[num_surfloc 1 3]);

%blend in the skin tone
skincol = .75*[0.9373 0.8157 0.8118];
C(:,1,1) = C(:,1,1).*blend_alpha + skincol(1)*(1-blend_alpha);
C(:,1,2) = C(:,1,2).*blend_alpha + skincol(2)*(1-blend_alpha);
C(:,1,3) = C(:,1,3).*blend_alpha + skincol(3)*(1-blend_alpha);

set(hsurf,'FaceVertexCData',squeeze(C));

caxis([surf_min surf_max]);
colormap(cmap);
%%%%% remove extra contour points, if around
if(~isempty(hcontour))
    for i = 1:numel(hcontour)
        loc = [get(hcontour(i),'XData'); get(hcontour(i),'YData'); get(hcontour(i),'ZData')]';
        loc_incap = inhull(loc,eeg_loc,khull);
        set(hcontour(i),'XData',loc(loc_incap,1));
        set(hcontour(i),'YData',loc(loc_incap,2));
        set(hcontour(i),'ZData',loc(loc_incap,3));
    end
end

%%%%%
if(~show_electrodes)
    delete(hline);
    delete(htext);
end

set(ax,'userdata',userdata);