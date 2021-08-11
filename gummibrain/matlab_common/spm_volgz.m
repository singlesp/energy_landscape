function [V vol] = spm_vol_gz(fname,varargin)

V=[];
vol = [];

inputchar=false;
if(ischar(fname))
    inputchar=true;
    fname={fname};
end


for f = 1:numel(fname)
    [~,~,ext]=fileparts(fname{f});
    if(strcmpi(ext,'.gz'))
        tmp=tempdir;
        ftmp=gunzip(fname{f},tmp);
        ftmp=ftmp{1};
        Vtmp = spm_vol(ftmp);
        voltmp=spm_read_vols(Vtmp,varargin{:});
        delete(ftmp);
    else
        Vtmp = spm_vol(fname{f});
        voltmp=spm_read_vols(Vtmp,varargin{:});
    end
    Vtmp.fname=fname{f};
    if(isempty(V))
        V=Vtmp;
        vol = {voltmp};
    else
        V(end+1) = Vtmp;
        vol{end+1}=voltmp;
    end
end

if(inputchar && numel(fname) == 1)
    vol = vol{1};
end
