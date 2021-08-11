function [ cifti ] = ciftiseparate(filename,caret7command)
%Open a CIFTI file by converting to GIFTI external binary first and then
%using the GIFTI toolbox

if(~exist('caret7command','var') || isempty(caret7command))
    [~,caret7command]=system('which wb_command');
    caret7command=regexp(caret7command,'\n','split');
    caret7command=caret7command(~cellfun(@isempty,caret7command));
    caret7command=caret7command{1};
end

grot=fileparts(filename);
if (size(grot,1)==0)
grot='.';
end
tmpname = tempname(grot);

tmpname_nii=[tmpname '.nii'];
tmpname_lh=[tmpname '_lh.gii'];
tmpname_rh=[tmpname '_rh.gii'];

unix(sprintf('%s -cifti-separate %s COLUMN -volume-all %s -metric CORTEX_LEFT %s -metric CORTEX_RIGHT %s',...
    caret7command, filename,tmpname_nii,tmpname_lh,tmpname_rh));

nii=read_avw(tmpname_nii);
lh_gii=gifti(tmpname_lh);
rh_gii=gifti(tmpname_rh);

cifti=struct('volume',nii,'lh',lh_gii.cdata,'rh',rh_gii.cdata);

unix(sprintf('rm %s %s %s',tmpname_nii,tmpname_lh,tmpname_rh));

end


