function c = set_combinations(sets)
%c = set_combinations(set_sizes)
%
% generates a collection of all combinations
% of elements in lists given by cell array sets.
% 
% c = [N x 3] where N = (set1 size)*(set2 size)*(set3 size)...
%
% ex:
% set_combinations({1, 10:11, 100:102})
% 
% ans =
% 
%      1    10   100
%      1    11   100
%      1    10   101
%      1    11   101
%      1    10   102
%      1    11   102

s = cellfun(@numel,sets);
n = numel(s);

c = zeros(prod(s),n);

for i = 1:n
    a = 1:prod(s(1:i));
    a = fix((a-1)/prod(s(1:i-1)))+1;
    idx = repmat(a,1,prod(s(i+1:end)));
    c(:,i) = sets{i}(idx)';
end
