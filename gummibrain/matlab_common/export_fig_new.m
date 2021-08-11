function [im alpha] = export_fig_new(filename,varargin)
if(nargin < 1)
    [im alpha] = export_fig();
else
    if(~isempty(filename))
        [im alpha] = export_fig(filename,varargin{:});
    end
end
