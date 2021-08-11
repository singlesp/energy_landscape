function niifiles = analyze2nii(imgfiles,varargin)

p = inputParser;

p.addParamValue('merge4d',false);
p.addParamValue('cleanup',false);
p.addParamValue('mergename',[]);

p.parse(varargin{:});
r = p.Results;

do_merge = r.merge4d;
do_cleanup = r.cleanup;
nii_mergename = r.mergename;

if(ischar(imgfiles))
    if(size(imgfiles,1) == 1)
        imgfiles = {imgfiles};
    else
        imgfiles = cellstr(imgfiles);
    end
end

isimg = regexpimatch(imgfiles,'\.img$');
hdrfiles = regexprep(imgfiles,'.\img$','.hdr');

niifiles = {};
if(do_merge)
    Vhdr = spm_vol(char(imgfiles));
    [d,fname,ext] = fileparts(imgfiles{1});
    if(~isempty(nii_mergename))
        [~,fname,~] = fileparts(nii_mergename);
    end
    niifile = fullfile(d,[fname '.nii']);
    spm_file_merge_kj(char(imgfiles),niifile);
    
    success = false;
    if(exist(niifile,'file'))
        Vnii = spm_vol(niifile);
        if(numel(Vnii) == numel(Vhdr))
            success = true;
        end
    end
    if(success)
        niifiles = {niifile};
        if(do_cleanup)
            delete(imgfiles{:});
            delete(hdrfiles{isimg});
        end
    end
else
    niifiles = {};
    for i = 1:numel(imgfiles)
        [d,fname,ext] = fileparts(imgfiles{i});
        niifile = fullfile(d,[fname '.nii']);
        Vhdr = spm_vol(imgfiles{i});
        V = spm_read_vols(Vhdr);
        Vhdr.fname = niifile;
        spm_write_vol(Vhdr,V);
        success = false;
        if(exist(niifile,'file'))
            Vnii = spm_vol(niifile);
            if(numel(Vnii) == numel(Vhdr))
                success = true;
            end
        end        

        if(success)
            niifiles{end+1} = niifile;
            if(do_cleanup)
                delete(imgfiles{i});
                if(isimg(i))
                    delete(hdrfiles{i});
                end
            end
        end
    end
end