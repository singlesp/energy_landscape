function labelinfo = mri_ribbonlabels(mrifile)
% labelinfo = mri_ribbonlabels(mrifile)
% 
% mrifile: can be structural filename, spm_vol struct, or actual 3d volume image
% labelinfo: struct containing label values gm_lh,wm_lh,gm_rh,wm_rh
%
% if mrifile='default', just return default ribbon label values
%
if(ischar(mrifile) && strcmpi(mrifile,'default'))
    gm_rh = 42;
    wm_rh = 41;
    gm_lh = 3;
    wm_lh = 2;
    labelinfo = fillstruct(gm_lh,wm_lh,gm_rh,wm_rh);
    return;
end

if(ischar(mrifile))
    Vhdr = spm_vol(mrifile);
    V = spm_read_vols(Vhdr);
elseif(isstruct(mrifile))
    V = spm_read_vols(mrifile);
elseif(isnumeric(+mrifile))
    V = +mrifile;
end


anatlabels = unique(V(:));
anatlabels = anatlabels(anatlabels ~= 0);

if(numel(anatlabels) > 4)
    warning('ribbon file must have exactly 4 labels. quitting.');
    labelinfo = [];
    return;
end
%%
label_count = zeros(size(anatlabels));
label_bounds = zeros(numel(anatlabels),2);
for i = 1:numel(anatlabels)
    Vi = V==anatlabels(i);
    label_count(i) = sum(Vi(:));
    zmip = squeeze(max(Vi,[],2));
    xzmip = max(zmip,[],2);
    yzmip = max(zmip,[],1);
    
    xzmip_idx = find(xzmip > 0);
    label_bounds(i,:) = [min(xzmip_idx) max(xzmip_idx)];
end

%%
[~,idx1] = sort(label_bounds(:,1),'ascend');
[~,idx2] = sort(label_bounds(:,2),'ascend');

gm_rh = anatlabels(idx1(1));
wm_rh = anatlabels(idx1(2));

gm_lh = anatlabels(idx2(end));
wm_lh = anatlabels(idx2(end-1));

labelinfo = fillstruct(gm_lh,wm_lh,gm_rh,wm_rh);
