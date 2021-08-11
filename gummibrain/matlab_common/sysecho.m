function v = sysecho(varargin)

fprintf('%s\n',varargin{1});
v = system(varargin{:});
