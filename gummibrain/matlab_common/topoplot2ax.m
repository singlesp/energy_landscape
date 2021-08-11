function varargout = topoplot2ax(ax,varargin)
axes(ax);
cla(ax);
varargout = topoplot(varargin{:});
axis tight;