function d = localtoolboxdir

if(isunix)
    d = '/home/kjamison/TOOLBOX';
elseif(ispc)
    %d = 'E:/Users/Keith/MATLAB_Toolboxes';
    d = 'D:/MATLAB_Toolboxes';
else
    warning('unknown toolbox platform');
end
