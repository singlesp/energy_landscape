%% plot centroids on gummibrain

clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);

%% Add paths for FSL, SPM, and keith's matlab_common

if(isempty(which('read_avw')))
    setenv('FSLDIR','/usr/local/fsl');
    addpath('/usr/local/fsl/etc/matlab');
end
if(isempty(which('spm_vol')))
    addpath('~/Documents/MATLAB/spm12');
end
%if(isempty(which('iso2meshver')))
%    addpath('~/MATLAB_TOOLBOXES/iso2mesh');
%end
if(isempty(which('clearb')))
    addpath('~/Documents/gummibrain/matlab_common');
end


%% plot f

numClusters = 4;
split=22; % i am using split to denote different processing applied to data 
% split 0 corresponds to gsr only (LSDgsr_cat.mat)
% split 1 corresponds to gsr + bp 0.008-0.09 (LSDgsr_bp_cat.mat)
% split 2 corresponds to gsr + bp 0.008-0.2 (LSDgsr_bp2_cat.mat)
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
nsubjs=15;
nscans=29;
TR=217;
[nparc,~] = size(centroids);

if nparc == 90
    atlas_file = 'aal_space-MNI152NLin6_res-2x2x2';
    
elseif nparc == 232
    atlas_file = 'Schaefer2018_232Parcels_7Networks_order_FSLMNI152_2mm';
    
elseif nparc == 116
    atlas_file = 'Schaefer2018_116Parcels_7Networks_order_FSLMNI152_2mm';
    
elseif nparc == 454
    atlas_file = 'Schaefer2018_454Parcels_7Networks_order_FSLMNI152_2mm';
    
elseif nparc == 463 || nparc == 462
    atlas_file = 'lausanne2008_scale4_MNI152_463parcels';
elseif nparc == 461
    atlas_file = 'lausanne2008_scale4_MNI152_463parcels';
    load /Users/sps253/Documents/ROI_maps/Lausanne_463_subnetworks.mat subnetworks_reorder
%     subnetworks_reorder = NaN(463,1);
%     subnetworks_reorder(subnetworks==3)=1; %
%     subnetworks_reorder(subnetworks==2)=2; %
%     subnetworks_reorder(subnetworks==5)=3;
%     subnetworks_reorder(subnetworks==4)=4;
%     subnetworks_reorder(subnetworks==7)=5; %
%     subnetworks_reorder(subnetworks==6)=6; %
%     subnetworks_reorder(subnetworks==1)=7; %
%     subnetworks_reorder(subnetworks==8)=8; %
    
    centroids = vertcat(centroids(1:13,:),[0 0 0 0],centroids(14:end,:),[0 0 0 0]);
    roimask = repelem(1,463);
    roimask([14 463])=0;
    roimask = roimask';
    roimask(subnetworks_reorder==8)=0;
elseif nparc == 333
    atlas_file = 'Parcels_MNI_333';
end

%%

atlas_fullpath = ['/Users/sps253/Documents/parcs/', atlas_file];

blobs = make_atlas_blobs([atlas_fullpath, '.nii'],'volumesmoothing','true');
%set colormap
cmap=plasma;

%choose atlas
%options for whichatlas: aal, cc200, cc400, ez, ho, tt, fs86
whichatlas={'333'};
clc;

%%
%set data you want to plot
% load 5HTvecs_ls463
% data = horzcat(zscore(mean5HT2A_ls463),zscore(mean5HT1A_ls463),zscore(mean5HT1B_ls463),zscore(mean5HTT_ls463));
% data=event_centroids;
data = network7labels;
%set min/max limits for plot
clim=[0 7];

cmap = lines(7);
% cmap(7,:) = cmap(7,:)*0.5;
% cmap(8,:) = cmap(7,:)*0.5;
% cmap = plasma;
% 
f=figure;


for i=1:numClusters
    
    subplot(1,numClusters,i);
    
    img=display_atlas_blobs(data(:,i),blobs,...
        'atlasname',whichatlas,...
        'render',true,...
        'backgroundimage',true,...
        'crop',true,...
        'colormap',cmap,...
        'clim', clim,...
        'backgroundcolor',[1 1 1]);%,...
%         'roimask',roimask); %last argmument optional, depends if you need to remove roi's
    
    
    imshow(img);
    c=colorbar('SouthOutside', 'fontsize', 16);
%     c.Label.String='Low-Amplitude Activity  High Amplitude Activity';
    set(gca,'colormap',cmap);
    caxis(clim);
%     title(clusterNames{i});
    
%     annotation('textbox',[.89 .33 .1 .2],'String','RH','EdgeColor','none','fontname','Arial','fontsize',16,'color','black')
%     annotation('textbox',[.07 .33 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','black')
%     annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none','fontname','Arial', 'fontsize',16,'color','black', 'horizontalalignment','center')
%     annotation('textbox',[.45 .1 .1 .04],'String','Medial','EdgeColor','none','fontname','Arial', 'fontsize',16,'color','black', 'horizontalalignment','center')
%     
end

% figure;
% title(clusterNames{i});
% 