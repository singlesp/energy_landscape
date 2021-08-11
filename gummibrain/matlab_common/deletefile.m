function deletefile(varargin)

w = warning('off','MATLAB:DELETE:FileNotFound');
delete(varargin{:});
warning(w);
