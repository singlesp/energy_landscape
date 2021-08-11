function varargout = export_fig_and_close(fig,varargin)
varargout = {export_fig(varargin{:},fig)};
close(fig);
