function varargout = addglmdenoisepath

origpath = path;
[d,~,~] = fileparts(which('GLMdenoisedata'));
addpath(fullfile(d,'utilities'));

if(nargout > 0)
    varargout = {origpath};
end