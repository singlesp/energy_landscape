function newstr = fixslash(str)

newstr = regexprep(str,'\','\\\\');