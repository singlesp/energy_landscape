clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 4;
split=33; % i am using split to denote different processing applied to data 
% split 0 corresponds to gsr only (LSDgsr_cat.mat)
% split 1 corresponds to gsr + bp 0.008-0.09 (LSDgsr_bp_cat.mat)
% split 2 corresponds to gsr + bp 0.008-0.2 (LSDgsr_bp2_cat.mat)
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
% TR = 217;


% addpaths;
% masterdir = fullfile(basedir,'results',name_root);
% 
% %%
% load(fullfile(masterdir,'clusterAssignments',['k',num2str(numClusters),name_root,'.mat']));
% overallPartition = clusterAssignments.(['k',num2str(numClusters)]).partition;
% centroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
 %clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
% savedir = fullfile(masterdir,'analyses','centroids');
% mkdir(savedir);

%%

% centroids=init_centroids;
% overallNames = {'5-HT2a','5-HT1a','5-HT1b','5-HT4','5-HTT'};
% overallNames = {'MS-1a','MS-2a','MS-2a','MS-2b'};
% clusterNames = NAME_CLUSTERS_ANGLE(centroids); 
overallNames = clusterNames;
% centroids = horzcat((mean5HT2A_sch454/max(mean5HT2A_sch454)),(mean5HT1A_sch454/max(mean5HT1A_sch454)),(mean5HT1B_sch454/max(mean5HT1B_sch454)),(mean5HT4_sch454/max(mean5HT4_sch454)),(mean5HTT_sch454/max(mean5HTT_sch454)));%centroids;
% centroids = double(centroids > 0.85);
% centroids = horzcat(mean5HT2A_sch454,mean5HT1A_sch454,mean5HT1B_sch454,mean5HT4_sch454,mean5HTT_sch454);%centroids;
% centroids = centroids/(max(max(centroids)));
% centroids = squeeze(mean(PLcentroids(:,:,:)));

[nparc,numClusters] = size(centroids);
% 
[~,~,~,net7angle] = NAME_CLUSTERS_ANGLE(centroids);
% 
YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

%% plot

clusterColors = GET_CLUSTER_COLORS(numClusters);
%YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;247 253 205; 218 166 86; 199 109 117] / 255;
% YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117] / 255;
YeoColors = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;];
YeoColors = [YeoColors;YeoColors];

f = figure;
imagesc(net7angle); ax = gca;
colormap('plasma')
set(ax,'xaxisLocation','top')
xticks(1:14); xticklabels(YeoNetNames);
xtickangle(90);
for K = 1:14
	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	YeoColors(K,:), ax.XTickLabel{K});
end
yticks(1:numClusters); COLOR_TICK_LABELS(false,true,numClusters);
set(ax,'FontSize',8);

%%

[~,~,net7angle_Up,net7angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids);
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
numNets = numel(YeoNetNames);

%% plot

%-7/27/20 SOMETHING SEEMS OFF WITH THIS FIGURE - DOES NOT SEEM TO MATCH FIG
%1 OR RADIAL PLOTS. LOOK AT SI LATER TO INVESTIGATE!!

% sim = [net7angle_Up;net7angle_Down];
% reOrder = [1 6 2 7 3 8 4 9 5 10];   % order 1+ 1- 2+ 2- etc.
% clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
% clusterColors(6:10) = clusterColors;
% clusterColors = hex2rgb(clusterColors);
% sim = sim(reOrder,:);
% clusterColors = clusterColors(reOrder,:);
% yticklabs = {'1+','2+','3+','4+','5+','1-','2-','3-','4-','5-'};
% f = figure;
% imagesc(sim); ax = gca;
% colormap('plasma')
% set(ax,'xaxisLocation','top')
% xticks(1:(numNets*2)); xticklabels(YeoNetNames);
% xtickangle(90);
% for K = 1:(numNets*2)
% 	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
% 	YeoColors(K,:), ax.XTickLabel{K});
% end
% yticks(1:numClusters*2); yticklabels(yticklabs(reOrder)');
% for K = 1:10 %(numClusters*2)
% 	ax.YTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
% 	clusterColors(K,:), ax.YTickLabel{K});
% end
% set(ax,'FontSize',8);
% set(ax,'TickLength',[0 0]);
% box off
% 
% for i = 2:2:8
%     line([0 7.5],[i + 0.5 i + 0.5],'LineWidth',1,...
%         'color',[0.5 0.5 0.5]);
% end
% [~,sysIdx] = max(sim,[],2);
% for i = 1:length(sysIdx)
%     rectangle('Position',[sysIdx(i)-0.5,i-0.5,1,1],'EdgeColor','k',...
%         'LineWidth',2);
% end
% 
% f.PaperUnits = 'inches';
% f.PaperSize = [2.5 2];
% f.PaperPosition = [0 0 2.5 2];
% saveas(f,fullfile(savedir,['Systems_k',num2str(numClusters),name_root,'.pdf']));


%% make radial plots

clusterColors = GET_CLUSTER_COLORS(numClusters);

clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{8} = '';

% net7angle_Up=abs(net7angle_Up-net7angle_UpPL);
% net7angle_Down=abs(net7angle_Down-net7angle_DownPL);

f=figure;
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net7angle_Up(K,:) net7angle_Up(K,1)],'k');
    polarplot(netAngle,[net7angle_Down(K,:) net7angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rlim([0.0 0.7]);
    rticks([0.2 0.4 0.8]); rticklabels({'','0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];
% saveas(f,fullfile(savedir,['SystemsRadial_k',num2str(numClusters),name_root,'.pdf']));
% 
% % for source data file
% save(fullfile(savedir,['Fig2b__YeoSystemAlignment_k',num2str(numClusters),'.mat']),'netAngle','net7angle_Up','net7angle_Down');
