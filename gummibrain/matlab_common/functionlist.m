function functionlist(functionlist, varargin)
for i = 1:numel(functionlist)
    functionlist{i}(varargin{:});
end
