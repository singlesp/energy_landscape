%% make vectors for all 5-HT data

% 5-HT data is from Beliveau et al 2017 Neuroscience: High resolution in
% vivo atlas of human serotonin system

%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/5-HT_atlas';
cd(basedir);
%%

maxROI = 232;

if maxROI == 90
    atlas_file = 'aal_space-MNI152NLin6_res-2x2x2';
elseif maxROI == 116
    atlas_file = 'Schaefer2018_116Parcels_7Networks_order_FSLMNI152_2mm';
elseif maxROI == 232
    atlas_file = 'Schaefer2018_232Parcels_7Networks_order_FSLMNI152_2mm';
elseif maxROI == 454
    atlas_file = 'Schaefer2018_454Parcels_7Networks_order_FSLMNI152_2mm';
elseif maxROI == 463
    atlas_file = 'lausanne2008_scale4_MNI152_463parcels';
end


atlas_fullpath = ['/Users/sps253/Documents/parcs/', atlas_file];
nii = load_untouch_nii([atlas_fullpath, '.nii']);
voxel_atlas = nii.img;

labellist = setxor(unique(voxel_atlas), 0);
labellist = labellist(1:maxROI);

mean5HT2A_sch232 = double(ts_extract('5-HT2A/5HT2a_mean_bmax_2mmMNI',[atlas_file,'.nii']));
mean5HT1A_sch232 = double(ts_extract('5-HT1A/5HT1a_mean_bmax_2mmMNI',[atlas_file,'.nii']));
mean5HT1B_sch232 = double(ts_extract('5-HT1B/5HT1b_mean_bmax_2mmMNI',[atlas_file,'.nii']));
mean5HT4_sch232 = double(ts_extract('5-HT4/5HT4_mean_bmax_2mmMNI',[atlas_file,'.nii']));
mean5HTT_sch232 = double(ts_extract('5-HTT/5HTT_mean_bmax_2mmMNI',[atlas_file,'.nii']));

mean5HT1A_sch232(isnan(mean5HT1A_sch232))=0;

%%
savedir = '/Users/sps253/Documents/brain_states-master/';
cd(savedir);
save(fullfile([savedir,'5HTvecs_sch232.mat']),'mean5HT1A_sch232','mean5HT1B_sch232','mean5HT2A_sch232','mean5HT4_sch232','mean5HTT_sch232');