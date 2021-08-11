function x = flattencell(c)
%isc = cellfun(@iscell,c);
%a = cellfun(@flattencell,c(isc),'uniformoutput',false)
c = cellfun(@(x)(rot90(flatten(x))),c,'uniformoutput',false);
x = flatten([c{:}]);