function idx = StringIndex(all_names, names)
% idx = StringIndex(all_names, names)
% all_names = list of strings
% names = Nx1 cell array where each item can be:
%       a) a single string to find in all_names
%       b) a {'str1','str2'} array to include the entire range from str1:str2
% 
% e.g.
%   allnames = {'a','b','c','d','e','f'};
%   idx = StringIndex(allnames,{'b','c','f'}) --> [2 3 6]
%   namerange = {'c','f'};
%   idx = StringIndex(allnames,{'a',namerange}) --> [1 3:6]

if(~iscell(names))
    ix = find(cellfun(@(x)(strcmpi(x,names)),all_names));
    if(isempty(ix))
        idx = [];
    else
        idx = ix;
    end
else
    idx = [];
    for i = 1:numel(names)
        if(ischar(names{i}))
            ix = find(cellfun(@(x)(strcmpi(x,names{i})),all_names));
            if(isempty(ix))
                %ix = -1;
                ix = [];
            end
            %idx = [idx ix(1)];
            idx = [idx ix];
        %elseif(iscell(names{i}))
        %    ix1 = find(cellfun(@(x)(strcmpi(x,names{i}{1})),all_names));
        %    ix2 = find(cellfun(@(x)(strcmpi(x,names{i}{2})),all_names));
        end
    end
end
