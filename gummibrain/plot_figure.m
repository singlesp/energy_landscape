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


%% plot figures
atlasblobs_list=load('./atlasblobs_saved.mat');
atlasblobs_list=atlasblobs_list.atlasblobs;

%set colormap
cmap=hot;

%choose atlas
%options for whichatlas: aal, cc200, cc400, ez, ho, tt, fs86
whichatlas={'fs86'}
clc;

%set data you want to plot
data=[];

%set min/max limits for plot
clim=[];

img=display_atlas_blobs(data,atlasblobs_list,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'clim', clim);

figure;
imshow(img);

%to show colorbar
%c=colorbar('SouthOutside', 'fontsize', 16);
%c.Label.String='';
%set(gca,'colormap',cmap);
%caxis([0.5 8.5]);
%title('TT', 'fontsize',18);

% to annotate RH, LH, medial, lateral sides fo the brain
%annotation('textbox',[.835 .38 .1 .2],'String','RH','EdgeColor','none', 'fontsize',16,'color','white')
%annotation('textbox',[.12 .38 .1 .2],'String','LH','EdgeColor','none', 'fontsize',16,'color','white')
%annotation('textbox',[.45 .9 .1 .04],'String','Lateral','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center','backgroundcolor','black')
%annotation('textbox',[.45 .21 .1 .04],'String','Medial','EdgeColor','none', 'fontsize',16,'color','white', 'horizontalalignment','center', 'backgroundcolor','black')

%save image
%imwrite(img,sprintf('~/Downloads/%s_t.png',whichatlas));
