function newpath = uiputfile_fullpath(FilterSpec,DialogTitle,DefaultName)
v={FilterSpec};
if(exist('DialogTitle','var'))
    v{end+1}=DialogTitle;
    if(exist('DefaultName','var'))
        v{end+1}=DefaultName;
    end
end
[f,p]=uiputfile(v{:});
if(isnumeric(f) && f==0)
    newpath='';
else
    newpath=fullfile(p,f);
end
