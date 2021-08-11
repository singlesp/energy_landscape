function [deplist filedates] = DependencyInfo(mfile,sortby,varargin)

if(nargin < 1 || isempty(mfile))
    db = dbstack(1,'-completenames');
    mfile = db.file;
end

if(nargin < 2 || isempty(sortby))
    sortby = 'date';
end

deplist = depfun(mfile,'-toponly','-quiet');

deplist = deplist(cellfun(@isempty,regexpi(deplist,fixslash(['^' matlabroot]))));
%deplist = deplist(~strcmp(which(mfile),deplist));
filedates = [];
for i = 1:numel(deplist)
    d=dir(deplist{i});
    filedates(i) = d.datenum;
end

if(strcmpi(sortby,'date'))
    [~,sortidx] = sort(filedates,'descend');
else
    [~,sortidx] = sort(deplist);
end


deplist = deplist(sortidx);
filedates = filedates(sortidx);

for i = 1:numel(deplist)
    fprintf('%-25s %s\n',datestr(filedates(i),'yyyy-mm-dd HH:MM:SS PM'), deplist{i});
end