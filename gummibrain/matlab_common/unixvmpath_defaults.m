function dirmaps = unixvmpath_defaults
%ex:
% dirmaps = {'C:/folder/folder2', '/cygdrive/c/folder/folder2'}

%VirtualBox also mounts the users Temp dir, so make any references
% to "c:\....\Temp" or "/tmp" to that mount point    
tmpd = regexprep(tempdir,'[/\\]+$','');    

dirmaps = {
    researchdir,                '/home/brain/HostData/Data', ...
    {tmpd, tempdir, '/tmp'} ,   '/home/brain/HostTemp', ...
    'D:/',                      '/home/brain/HostData/'
};

