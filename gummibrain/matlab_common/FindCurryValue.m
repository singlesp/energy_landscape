function val = FindCurryValue(fieldnames, values, whichname)
% val = FindCurryValue(fieldnames, values, whichname)

val = values{find(strcmpi(fieldnames,whichname))};
