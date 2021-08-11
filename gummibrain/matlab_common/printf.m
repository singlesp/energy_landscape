function varargout = printf(varargin)
output = fprintf([varargin{1} '\n'],varargin{2:end});
if(nargout > 0)
    varargout = {output};
end