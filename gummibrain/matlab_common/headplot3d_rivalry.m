function ax = headplot3d(subj,vals,which_head,show_labels,do_skinblend,splinevals,minmax,existing_ax,cmap)
%
%
% ax = headplot3d(subj,vals,which_head,show_labels,do_skinblend,splinevals,minmax)
%
%
% subj = subject name
% vals = electrode values (or blank if providing precalculated surface
%   interpolation (splinevals)
%
% which_head = ('default','skin','sphere') ... which head surface to use
%   [default = 'skin']
%
% show_labels = show electrode names? [default = false]
%
% do_skinblend = paint parts of the head outside the cap skin-colored?
%   (instead of interpolating values on the entire head) [default = true]
%
% splinevals = pre-calculated spline surface values (instead of
%   interpolating as you display)
%
% minmax = [minval maxval] for color scaling

if(nargin < 3 || isempty(which_head))
    which_head = 'skin';
end
which_head = lower(which_head);

if(nargin < 4)
    show_labels = false;
end

if(nargin < 5 || do_skinblend == true)
    blend = 'blend';
else
    blend = 'noblend';
end

if(nargin < 6)
    splinevals = [];
end

if(nargin < 7)
    minmax = [];
end

if(nargin < 8)
    existing_ax = [];
end

if(nargin < 9)
    cmap = [];
end

[bad_elec good_elec] = GetBadElectrodes128(subj);

%do spline interpolation if not already provided
if(isempty(splinevals))
    splinefile = GetSplineFilename(subj,numel(vals),which_head);
    splinevals = EvaluateSplineSurface(vals,splinefile);
end

reset = isempty(existing_ax);

if(reset)
    %first plot a blank head to fill in later
    %headplot converts values in RGB colors, which makes it difficult to
    %rescale, so let it plot blank values, then fill in with numerical values
    %afterward
    switch which_head
        case 'default'
            splinefile = sprintf('%s/%s_128/eeglab_spline_%d.mat',rivalrydir,subj,numel(vals));
            ax = headplot(zeros(size(good_elec)),splinefile,'labels',2);
        case {'skin','subject'}
            splinefile = sprintf('%s/%s_128/%s_spline_%d.mat',rivalrydir,subj,subj,numel(vals));
            skinfile = sprintf('%s/%s_128/%s_skin.mat',rivalrydir,subj,subj);
            ax = headplot(zeros(size(vals)),splinefile,'meshfile',skinfile,'labels',2);
        case 'sphere'
            splinefile = sprintf('%s/%s_128/sphere_spline_%d.mat',rivalrydir,subj,numel(vals));
            skinfile = sprintf('%s/%s_128/%s_sphere.mat',rivalrydir,subj,subj);
            ax = headplot(zeros(size(good_elec)),splinefile,'meshfile',skinfile,'labels',2);
    end
else
    ax = existing_ax;
end

%fill in surface with numerical values
[surf_vals valid_idx] = GetVertexColors(ax);
new_surf_vals = zeros(size(surf_vals,1),1);
new_surf_vals(valid_idx) = splinevals;
SetVertexColors(ax,new_surf_vals);

clean_headplot(ax,show_labels,blend,minmax,reset,cmap);

if(reset)
    axis vis3d equal;
    view([0 60]);
end


%topoplot(vals,eloc,'electrodes','labelpoint');