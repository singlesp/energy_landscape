%% Step 0: Add paths for FSL, SPM, and keith's matlab_common

if(isempty(which('read_avw')))
    setenv('FSLDIR','/usr/local/fsl');
    addpath('/usr/local/fsl/etc/matlab');
end
if(isempty(which('spm_vol')))
    addpath('~/MATLAB_TOOLBOXES/spm12');
end
%if(isempty(which('iso2meshver')))
%    addpath('~/MATLAB_TOOLBOXES/iso2mesh');
%end
if(isempty(which('clearb')))
    addpath('~/Source/matlab_common');
end

%% Step 1: build atlasblob structure for each atlas using make_atlas_blobs()
% each atlas might need tweaking in terms of extraction and smoothing 
% after building, an "atlasblobs" struct for each atlas, save this list so
% future calls to display_atlas_blobs() can just load existing struct and
% go

atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs;

%fsaverage:
%fs_default86.txt = created from merging FreeSurferROIlist86_Nypipe.txt with ColorLUT.txt
%mri_convert ${FREESURFER_HOME}/subjects/fsaverage/mri/aparc+aseg.mgz fsaverage_aparc+aseg.nii.gz
%flirt -in fsaverage_aparc+aseg.nii.gz -ref fsaverage_aparc+aseg.nii.gz -out fsaverage_aparc+aseg_2mm.nii.gz -applyisoxfm 2 -interp nearestneighbour
%labelconvert fsaverage_aparc+aseg_2mm.nii.gz ${FREESURFER_HOME}/FreeSurferColorLUT.txt fs_default86.txt fsaverage_nodes86_2mm.nii.gz

%load mid-sagittal slice from this file for displaying in background
bgfile=sprintf('%s/data/standard/MNI152_T1_1mm.nii.gz',getenv('FSLDIR'));


atlasdir='~/Downloads/atlases';


%whichatlas_all={'cc400','cc200','aal','ez','ho','tt','fs86'};
whichatlas_all={'cc400'};

atlasblobs_all=[];

for w = 1:numel(whichatlas_all)
    whichatlas=whichatlas_all{w};
  
    %smoothing=0;
    %smoothing=5;
    %smoothing=1;
    %volsmooth=false;
    
    %volsmooth = true/false do we perform a gentle voxel smoothing on each 
    % ROI before extracting surface? (helpful especially for small ROIs
    
    %surfsmoothing = iterations of surface vertex smoothing after
    %extracting from volume. Less necessary if volume smoothing performed
    
    atlasfile=sprintf('%s/%s_featimp.nii.gz',atlasdir,whichatlas);
    
    switch(lower(whichatlas))
        case 'cc400'
            surfsmoothing=0;
            volsmooth=true;
        case 'cc200'
            surfsmoothing=0;
            volsmooth=true;
        %case 'aal'
            %atlasfile=sprintf('%s/aal_roi_atlas116.nii.gz',atlasdir);
            %volsmooth=true;
            %surfsmoothing=1;
        %case 'fs86'
            %atlasfile='~/Downloads/fsaverage_nodes86_2mm.nii.gz';
            %volsmooth=true;
            %surfsmoothing=0;
        otherwise
            surfsmoothing=1;
            volsmooth=true;
            
    end
    
    atlasblobs=make_atlas_blobs(atlasfile,'backgroundfile',bgfile,...
        'volumesmoothing',volsmooth,'surfacesmoothing',surfsmoothing);

    atlasblobs.atlasname=whichatlas;
    
    if(isempty(atlasblobs_all))
        atlasblobs_all=atlasblobs;
    else
        atlasblobs_all(w)=atlasblobs;
    end
    
end

atlasblobs=atlasblobs_all;
save('~/Downloads/atlasblobs_saved.mat','atlasblobs');


%% Step 2: Use saved atlasblob surfaces to quickly render a given atlas

% for T stat 

atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs;


%whichatlas_all={'cc400','cc200','aal','ez','ho','tt','fs86'};
whichatlas_all={'fs86'}
close all;
clc;

for w = 1:numel(whichatlas_all)
    whichatlas=whichatlas_all{w};
    
    
    %for testing just assign y-axis position (AP axis) as value for each ROI
    %roivals=atlasblobs_list(strcmpi({atlasblobs_list.atlasname},whichatlas)).roicenters(:,2);
    %roivals=rand(size(roivals));
    %cmap=colormap_videen(256);
    %cmap=flipud(hot(256));
    %cmap=cmap(129:256,:);
    %cmap=colormap_hotcold(2000);
    %for i=786:1214
        %cmap(i,:)=[0.5, 0.5, 0.5];
    %end
    cmap=copper(256);
    
    img=display_atlas_blobs(temp',atlasblobs_list,...
        'atlasname',whichatlas,...
        'render',true,...
        'backgroundimage',false,...
        'crop',true,...
        'colormap',cmap,...
        'clim', [-6, 6]);
    %rescale(mean_fi(2,:).^2)
    %'alpha',rescale(pinv(rescale(abs(cc400_fi')))).^2
    %'alpha',cc400_t_facealpha
    figure;
    imshow(img);
    c=colorbar('SouthOutside', 'fontsize', 16);
    c.Label.String='TT: t-Statistic (Males > Females)';
    set(gca,'colormap',cmap);
    caxis([-10, 10]);
    %title('TT', 'fontsize',18); 
    %annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',16,'color','white')
    %annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','white')
    %annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center','backgroundcolor','black')
    %annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')
    
end
%imwrite(img,sprintf('~/Downloads/%s_t.png',whichatlas));


%% for feature importance

atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs;


%whichatlas_all={'cc400','cc200','aal','ez','ho','tt','fs86'};
whichatlas_all={'tt'}
close all;
clc;

for w = 1:numel(whichatlas_all)
    whichatlas=whichatlas_all{w};
    
    
    %for testing just assign y-axis position (AP axis) as value for each ROI
    %roivals=atlasblobs_list(strcmpi({atlasblobs_list.atlasname},whichatlas)).roicenters(:,2);
    %roivals=rand(size(roivals));
    %cmap=colormap_videen(256);
    %cmap=flipud(hot(256));
    %cmap=cmap(129:256,:);
    cmap=hot(256);
    cmap=cmap(1:254,:);
    
    img=display_atlas_blobs(rescale(abs(tt_fi')),atlasblobs_list,...
        'atlasname',whichatlas,...
        'render',true,...
        'backgroundimage',true,...
        'crop',true,...
        'colormap',cmap,...
        'clim', [0 1],...
        'alpha',rescale(pinv(rescale(abs(tt_fi')))));
    %rescale(mean_fi(2,:).^2)
    %
    %'alpha',cc400_t_facealpha
    figure;
    imshow(img);
    c=colorbar('SouthOutside', 'fontsize', 16);
    c.Label.String='TT: Feature Importance';
    set(gca,'colormap',cmap);
    caxis([0 1]);
    %title('TT', 'fontsize',18); 
    annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',16,'color','white')
    annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','white')
    annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center','backgroundcolor','black')
    annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')
    
end
%imwrite(img,sprintf('~/Downloads/%s_t.png',whichatlas));

%% for correlation between ROI size and HE

atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs;

%for testing just assign y-axis position (AP axis) as value for each ROI
%roivals=atlasblobs_list(strcmpi({atlasblobs_list.atlasname},whichatlas)).roicenters(:,2);
%roivals=rand(size(roivals));
%cmap=colormap_videen(256);
cmap=jet;


whichatlas={'fs86'}
clc;

img=display_atlas_blobs(temp,atlasblobs_list,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'clim', [-4 4]);

figure;
imshow(img);
%c=colorbar('SouthOutside', 'fontsize', 16);
%c.Label.String='';
%set(gca,'colormap',cmap);
%caxis([0.5 8.5]);
%title('TT', 'fontsize',18);
annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',16,'color','white')
annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','white')
annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center','backgroundcolor','black')
annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')


%imwrite(img,sprintf('~/Downloads/%s_t.png',whichatlas));
