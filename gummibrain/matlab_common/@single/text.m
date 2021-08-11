function h = text(varargin)
%w = warning('off','MATLAB:dispatcher:nameConflict');

isval = cellfun(@isnumeric,varargin);
varargin(isval) = cellfun(@double,varargin(isval),'uniformoutput',false);

h = builtin('text',varargin{:});

%warning(w.state,'MATLAB:dispatcher:nameConflict');
