function savehistory

histfile = sprintf('%s/history.m',prefdir);
if(~exist(histfile,'file'))
    histfile = sprintf('%s/History.xml',prefdir);
    if(~exist(histfile,'file'))
        return
    end
end

timestr = datestr(now,'yyyymmdd-HHMMss');
compstr = computerid;

histdir = [researchdir '/matlabhistory'];
if(~exist(histdir,'dir'))
    mkdir(histdir);
end

newfile = sprintf('%s/history_%s_%s.m',histdir,timestr,compstr);
copyfile(histfile,newfile,'f');
