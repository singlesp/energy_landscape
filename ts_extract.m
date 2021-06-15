function [ts] = ts_extract(image_file,atlas_file)
% NOTE: see comments on line 34 if not extracting enough regions. Also
% atlas file reads best if it is pasted in the baseDIR you are running from
% original function written by Elvisha Dhamala, PhD

%compute ts from a scan

V=read_avw(sprintf('%s',image_file));

%reshape to voxel x time series
V1=reshape(V,[],size(V,4));

%load in atlas file
atlas=read_avw(sprintf('%s',atlas_file));
atlas=reshape(atlas,[],1);


%combine atlas data
Vparc(:,1)=atlas;

%remove all voxels not present in any of the atlas masks
zero_atlas = find(all(Vparc==0,2));
V_nonzero=V1;
V_nonzero(zero_atlas,:)=[];
Vparc(zero_atlas,:)=[];

%get rid of all zero voxels still present in the time series data
zero_ts = find(all(V_nonzero==0,2));
V_nonzero(zero_ts,:)=[];
roi=unique(Vparc(:,1));
Vparc(zero_ts,:)=[];



%get ROIs from atlas
% roi=unique(Vparc(:,1)); % comment out this line and move it before line
% 29/30 if you are extracting less regions than you intended. If the result
% is the regions are extracted as NaN then your .nii file does not contain
% those regions
roi=roi(roi~=0);

for j=1:length(roi)
    %get time series for roi
    Vparc_roi=Vparc(:,1)==roi(j);
    V_roi=(V_nonzero(Vparc_roi,:));
    roi_ts(j,:)= nanmean(V_roi); %changed to nanmean so that regions with some nan voxels still get extracted
    
    
    clear Vparc_roi V_roi;
    
end
    
ts=roi_ts;

end

