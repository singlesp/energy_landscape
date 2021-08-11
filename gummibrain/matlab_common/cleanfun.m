function func_clean = cleanfun(func)
f = functions(func);
ws = f.workspace{1};
dumpstruct(ws)
fstrname_cleanfun = f.function;
clear f ws func;

func_clean = eval(fstrname_cleanfun);
clear fstrname_cleanfun;
