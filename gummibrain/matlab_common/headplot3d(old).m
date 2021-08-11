function ax = headplot3d(subj,vals,which_head,show_labels,do_skinblend,newpath)
% ax = headplot3d(subj,vals,which_head,show_labels,do_skinblend) 

if(nargin < 3)
    which_head = 'default';
end
which_head = lower(which_head);

if(nargin < 4)
    show_labels = false;
end

if(nargin < 5)
    do_skinblend = true;
end

if(nargin < 6)
    newpath = '';
end

if(show_labels)
    labels = 2;
else
    labels = 0;
end

% %plotting to subject head
% if(use_default_head)
%     splinefile = sprintf('%s_128/eeglab_spline_%d.mat',subj,numel(vals));
%     ax = headplot(vals,splinefile,'labels',labels);
% else
%     splinefile = sprintf('%s_128/%s_spline_%d.mat',subj,subj,numel(vals));
%     skinfile = sprintf('%s_128/%s_skin.mat',subj,subj);
%     ax = headplot(vals,splinefile,'meshfile',skinfile,'labels',labels);
% end

switch which_head
    case 'default'
        splinefile = sprintf('%s%s_128/eeglab_spline_%d.mat',newpath,subj,numel(vals));
        ax = headplot(vals,splinefile,'labels',labels);
    case {'skin','subject'}
        splinefile = sprintf('%s%s_128/%s_spline_%d.mat',newpath,subj,subj,numel(vals));
        skinfile = sprintf('%s%s_128/%s_skin.mat',newpath,subj,subj);
        ax = headplot(vals,splinefile,'meshfile',skinfile,'labels',labels);
    case 'sphere'
        splinefile = sprintf('%s%s_128/sphere_spline_%d.mat',newpath,subj,numel(vals));
        skinfile = sprintf('%s%s_128/%s_sphere.mat',newpath,subj,subj);
        ax = headplot(vals,splinefile,'meshfile',skinfile,'labels',labels);
end

ch = get(ax,'children');
hsurf = ch(strcmpi(get(ch,'type'),'patch'));

surf_val = get(hsurf,'FaceVertexCData');

%rescale the values to fit the original data range
sv_range = max(surf_val)-min(surf_val);
v_range = max(vals)-min(vals);
new_val = (v_range/sv_range)*(surf_val - min(surf_val)) + min(vals);
set(hsurf,'FaceVertexCData',new_val);

axis vis3d equal;
view([0 70]);

if(do_skinblend)
    clean_headplot(ax,true);
end

%topoplot(vals,eloc,'electrodes','labelpoint');