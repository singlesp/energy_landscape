function varargout = swap(a,b)
if(nargout == 0)
    if(isempty(inputname(1)) && ischar(a))
        varname1 = a;
        varname2 = b;
    else
        varname1 = inputname(1);
        varname2 = inputname(2);
    end
    
    tmp = sprintf('%0.5f',rand(1));
    tmpvarname = ['tmpvar' datestr(now,'HHMMSSFFF') tmp(3:end)];
    
    evalin('caller',[tmpvarname ' = ' varname1 ';']);
    evalin('caller',[varname1 ' = ' varname2 ';']);
    evalin('caller',[varname2 ' = ' tmpvarname ';']);
    evalin('caller',['clear ' tmpvarname ';']);
else
    varargout = {b,a};
end
