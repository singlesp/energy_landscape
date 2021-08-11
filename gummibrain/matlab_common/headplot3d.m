function ax = headplot3d(vals,skinfile,splinefile,show_labels,do_skinblend,splinevals,minmax,existing_ax,cmap,varargin)

% ax = headplot3d(vals,skinfile,splinefile,show_labels,do_skinblend,splinevals,minmax,existing_ax,cmap)
%
% vals = electrode values (or blank if providing precalculated surface
%   interpolation (splinevals)
%
% show_labels = show electrode names? [default = false]
%
% do_skinblend = paint parts of the head outside the cap skin-colored?
%   (instead of interpolating values over the entire head surface) [default = true]
%
% splinevals = pre-calculated spline surface values (instead of
%   interpolating as you display)
%
% minmax = [minval maxval] for color scaling.  
%          [val] will scale from [-abs(val) abs(val)]
%          default = [-abs(max(vals)) abs(max(vals))]
%
% existing_ax = if not empty, apply new vals to existing headplot.  This
%   can significantly speed up the drawing of a series of topographies
%   onto the same head model.
%
% cmap = colormap to be used for the headplot. eg: jet(100)
%
% %%%%%%%%%%%%%%%%%%%%%%
% Returns:
%   ax = axis the head was plotted to (can be passed to existing_ax next time)

chanind = [];
if(~isempty(varargin))
    idx = StringIndex(varargin,'chanind');
    if(~isempty(idx))
        chanind = varargin{idx+1};
        varargin = varargin([1:idx-1 idx+2:end]);
    end
end

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

if(nargin < 10)
    eloc = [];
end

%do spline interpolation if not already provided
if(isempty(splinevals))
    [splinevals Smat] = EvaluateSplineSurface(vals,splinefile);
    if(~isempty(chanind))
        vals = vals(chanind);
        Smat = Smat(:,chanind);
        splinevals = EvaluateSplineSurface(vals,Smat);
    end
end

%splinefile
%size(splinevals)

reset = isempty(existing_ax);

if(reset)
    
    if(isempty(vals))
        M = load(splinefile,'indices');
        vals = zeros(size(M.indices));
    end
    if(numel(vals) == 1)
        vals = zeros(vals,1);
    end
    
    %first plot a blank head to fill in later
    %headplot converts values in RGB colors, which makes it difficult to
    %rescale, so let it plot blank values, then fill in with numerical values
    %afterward
    if(isempty(skinfile))
        ax = headplot(zeros(size(vals)),splinefile,'labels',2);
    else
        ax = headplot(zeros(size(vals)),splinefile,'meshfile',skinfile,'labels',2);
    end
    hsurf = findobj(ax,'type','patch');
    set(hsurf,'tag','skin');
else
    ax = existing_ax;
    hsurf = findobj(ax,'type','patch','tag','skin');
end



%fill in surface with numerical values
[surf_vals valid_idx] = GetVertexColors(hsurf);

new_surf_vals = zeros(size(surf_vals,1),1);

new_surf_vals(valid_idx) = splinevals;
SetVertexColors(hsurf,new_surf_vals);


hsurf = clean_headplot(ax,show_labels,blend,minmax,reset,cmap,chanind);

if(~isempty(varargin))
	set(hsurf,varargin{:});
end

if(reset)
    axis vis3d equal;
    view([0 60]);
end
