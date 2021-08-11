function ROIViewer(xyz_mm, roi_info, varargin)

p = inputParser;
p.addParamValue('cvalues',[]);
p.addParamValue('clim',[]);
p.addParamValue('smooth',false);
p.addParamValue('parent',[]);
p.addParamValue('glassbrain',true);
p.addParamValue('roialpha',1);
p.addParamValue('showactive',false);
p.addParamValue('colormap',[]);
p.addParamValue('sort',false);
p.addParamValue('sortabs',false);

p.parse(varargin{:});
r = p.Results;

roi_cval = r.cvalues;
cl = r.clim;
do_smooth = r.smooth;
hfig = r.parent;
do_glassbrain = r.glassbrain;
roi_alpha = r.roialpha;
do_showactive = r.showactive;
cmap = r.colormap;
do_sort = r.sort;
do_sortabs = r.sortabs;

sort_type = 'descend';
if(ischar(do_sort))
    sort_type = do_sort;
    do_sort = true;
end

if(ischar(do_sortabs))
    sort_type = do_sortabs;
    do_sortabs = true;
end

if(do_showactive && ~isfield(roi_info,'activevoxels'))
    do_showactive = false;
end

if(isempty(cmap))
    cmap = jet(100);
end

if(isempty(roi_cval))
    roi_cval = [roi_info.cval];
end

if(isempty(cl))
    cl = [min(roi_cval) max(roi_cval)];
end

if(isempty(hfig))
    hfig = figure('color',[1 1 1]);
end


G = load('mni_glass_brain.mat');

lbroi = uicontrol(gcf,'style','listbox','units','normalized','position',[0 0 .25 1],...
    'min',1,'max',3,'callback',@roiList_click);

panright = uipanel('parent',gcf,'units','normalized','position',[.25 0 .75 1]);
axes('parent',panright);
hsrf = trisurf(G.Faces, G.Vertices(:,1), G.Vertices(:,2), G.Vertices(:,3),'facecolor',[0 0 1],'facealpha',.1,'linestyle','none','hittest','off');
axis vis3d equal;
camlight headlight;
lighting phong;
%shading interp;
RotationHeadlight;
%RotationKeypress;
hold on;

if(do_glassbrain)
    set(hsrf,'handlevisibility','off')
else
    delete(hsrf);
end


cval = roi_cval;
cval = (cval - min(cl))/(max(cl)-min(cl));


cval = max(1,min(size(cmap,1),round(cval*size(cmap,1))));
hp = [];
roi_descrip = [];
for i = 1:numel(roi_info)
    if(do_showactive)
        tmpxyz = xyz_mm(:,roi_info(i).activevoxels);
    else
        tmpxyz = xyz_mm(:,roi_info(i).voxels);
    end
    [roivert roitri] = roi2surf(tmpxyz,do_smooth);

    if(isempty(roitri) || do_showactive)
        hp(i) = plot3(tmpxyz(1,:),tmpxyz(2,:),tmpxyz(3,:),'.','color',cmap(cval(i),:));
    else
        hp(i) = patch('faces',roitri,'vertices',roivert,'linestyle','none','facecolor',cmap(cval(i),:),'facealpha',roi_alpha);
    end
    roi_descrip{i} = sprintf('(%d) %s : %.3f',roi_info(i).label, roi_info(i).name, roi_cval(i));
end

if(do_sort)
    [~,roi_sortidx] = sort(roi_cval,sort_type);
elseif(do_sortabs)
    [~,roi_sortidx] = sort(abs(roi_cval),sort_type);
else
    roi_sortidx = 1:numel(roi_info);
end
roi_descrip = roi_descrip(roi_sortidx);

set(hfig,'toolbar','figure');
lighting phong;
colormap(cmap);
colorbar;
set(gca,'clim',cl);
%legend(cleantext(roi_names));
%set(lbroi,'string',{roi_info.name});
set(lbroi,'string',roi_descrip);
setappdata(gcf,'appdata',fillstruct(roi_info,hp,lbroi,roi_sortidx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiList_click(hObj,event)
A = getappdata(gcbf,'appdata');
dumpstruct(A);
selidx = get(lbroi,'Value');
if(isempty(selidx))
    return;
end
selidx = roi_sortidx(selidx);

for i = 1:numel(roi_info)
    if(any(selidx == i))
        set(hp(i),'visible','on','hittest','on');
    else
        set(hp(i),'visible','off','hittest','off');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [verts faces] = roi2surf(roi_xyz,do_smooth)
v = roi_xyz';
if(size(v,1) < 4)
    verts = v;
    faces = [];
    return;
end

T = delaunayn(v);
[L ed] = edgelengths(T,v);
L = reshape(L,[],6);
longedge = any(L > 1.5*median(L(:)), 2) | any(L < .2*median(L(:)),2);
if(~all(longedge))
    T(longedge,:) = [];
end

T = surfedge(T);
T=meshreorient(v,T);


if(do_smooth)
    conn = meshconn(T,size(v,1));
    smoothiter = 3;
    useralpha = .5;
    usermethod = 'laplacian';
    userbeta = .5;
    v = smoothsurf(v,[],conn,smoothiter,useralpha,usermethod,userbeta);
end
%[v T] = meshcheckrepair(v,T,'deep');
% try
%     vmin = min(v,[],1);
%     v= v-repmat(vmin,size(v,1),1);
%     [v,T]=remeshsurf(v,[T ones(size(T,1),1)],median(L(:)));
%     v = v+repmat(vmin,size(v,1),1);
%     verts = v;
%     faces = T(:,1:3);
% catch
% %     
% end
verts = v;
faces = T;