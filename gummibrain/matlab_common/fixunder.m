function newstr = fixunder(str)

newstr = regexprep(str,'_','\\_');