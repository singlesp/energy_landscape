function rmglmdenoisepath(varargin)

if(nargin > 0)
    origpath = varargin{1};
    path(origpath);
else
    rmpath([justdir(which('GLMdenoisedata.m')) '/utilities']);
end