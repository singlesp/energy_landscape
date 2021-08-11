function v = SelectOutput(func,inputs,whichoutput)
allv = cell(1,max(whichoutput));
eval(['[ ' sprintf('allv{%d} ',1:max(whichoutput)) ' ] = func(inputs{:});']);
v = allv(whichoutput);
if(numel(whichoutput) == 1)
    v = v{1};
end
