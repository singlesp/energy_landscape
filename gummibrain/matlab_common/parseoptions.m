function options = parseoptions(defaultoptions,varargin)
% S = parseoptions(defaults, newoptions...)
%
% Inputs: 
%   defaultoptions: struct containing default parameters
%   newoptions: struct OR 'param',val list with options to overwrite
%               defaults.  If newoptions val=[], use default!!
%
% Example:
%  >> defaults = struct('firstparam',7,'secondparam',14,'thirdparam',23);
%  >> options = parseoptions(defaults,'firstparam',8,'thirdparam',[])
%  options = 
%     firstparam: 8
%     secondparam: 14
%     thirdparam: 23

if(isempty(defaultoptions))
    options=struct();
else
    options=defaultoptions;
end
newoptions=mergestruct(varargin{:});

fn=fieldnames(newoptions);
for f = 1:numel(fn)
    opt=newoptions.(fn{f});
    if(~(isnumeric(opt) && isempty(opt)))
        options.(fn{f})=newoptions.(fn{f});
    end
end
