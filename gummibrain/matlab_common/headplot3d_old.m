function ax = headplot3d(subj,vals,which_head,show_labels,do_skinblend)
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
        splinefile = sprintf('%s_128/eeglab_spline_%d.mat',subj,numel(vals));
        ax = headplot(vals,splinefile,'labels',labels);
    case {'skin','subject'}
        splinefile = sprintf('%s_128/%s_spline_%d.mat',subj,subj,numel(vals));
        skinfile = sprintf('%s_128/%s_skin.mat',subj,subj);
        ax = headplot(vals,splinefile,'meshfile',skinfile,'labels',labels);
    case 'sphere'
        splinefile = sprintf('%s_128/sphere_spline_%d.mat',subj,numel(vals));
        skinfile = sprintf('%s_128/%s_sphere.mat',subj,subj);
        ax = headplot(vals,splinefile,'meshfile',skinfile,'labels',labels);
end


axis vis3d equal;
view([0 60]);

if(do_skinblend)
    clean_headplot(ax,true);
end

%topoplot(vals,eloc,'electrodes','labelpoint');