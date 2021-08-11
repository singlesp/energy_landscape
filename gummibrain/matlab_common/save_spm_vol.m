function fname = save_spm_vol(data,hdr_template,varargin)
Vnew = struct(varargin{:});
if(~isfield(Vnew,'fname') || ~isempty(Vnew.fname))
    Vnew.fname = sprintf('%s.nii',tempname);
end
Vsave = hdr_template;
Vsave.descrip = '';
if(isfield(Vnew,'pinfo'))
    Vsave = rmfield(Vsave,'pinfo');
end

Vsave = mergestruct(Vsave,Vnew);
spm_write_vol(Vsave,data);

fname = Vsave.fname;
