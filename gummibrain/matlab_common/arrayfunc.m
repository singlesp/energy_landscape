function varargout = arrayfunc(varargin)
% [A,...] = arrayfunc(@func, array1,...)
%
% Handle arrayfunc intelligently.  Does not require 'uniformoutput', flag
% If outputs are all scalars return a matrix.  Otherwise return a cell array.
% 
allv = cell(1,nargout);
evstr = ['[ ' sprintf('allv{%d} ',1:nargout) ' ] = arrayfun(varargin{:},''uniformoutput'',false);'];
eval(evstr);

for i = 1:numel(allv)
    if(all(cellfun(@isscalar,allv{i})))
        allv{i} = cell2mat(allv{i});
%     elseif(all(cellfun(@isnumeric,allv{i})))
%         sz = size(allv{i}{1});
%         if(all(cellfun(@(x)(ndims(x)==numel(sz) && all(size(x)==sz)),allv{i})))
%             allv{i} = cell2mat(allv{i});
%        end
    end
    varargout{i} = allv{i};
end
