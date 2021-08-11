function r = compare_struct(s1,s2)

r = true;

if(isstruct(s1))
    fn = fieldnames(s1);
    for i = 1:numel(fn)
        disp(['***' fn{i}]);
        compare_struct(s1.(fn{i}),s2.(fn{i}));
    end
else
    if(s1 ~= s2)
        disp(s1);
        disp(s2);
    end
end