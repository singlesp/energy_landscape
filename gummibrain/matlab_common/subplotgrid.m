function h = subplotgrid(N,n,varargin)

c = ceil(sqrt(N));
r = ceil(N/c);

h = subplot(r,c,n,varargin{:});