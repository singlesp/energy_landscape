%% make figure 4a/b, i

clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);

%% set inputs
numClusters = 4;
split='main'
load(fullfile(['data/',split,'.mat']))
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

%% transition probabilities

% transition probability matrices:
% restTransitionProbabilityMats: each element reflects probability that
% state j occurs at the TR after state i, given that you are in state i (and not in the last TR of a scan)
% restTransitionProbabilityMatsNoPersist: each element reflectts
% probability that state j follows state i in the next state change,
% given that you are in state i (and not in the last TR of a scan). This
% *NoPersist version is used in the manuscript but both are interesting.

% get resting state transition probability matrices -- 2D is flattened by row for
% regressions


LSD_stop=TR*nscans;
PL_start=LSD_stop+1;
PL_stop=LSD_stop*2;


[LSDTransitionProbability2D,LSDTransitionProbabilityMats] = GET_TRANS_PROBS(partition(1:LSD_stop),subjInd(1:LSD_stop));
[LSDTransitionProbabilityNoPersist2D,LSDTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(partition(1:LSD_stop),subjInd(1:LSD_stop));

[PLTransitionProbability2D,PLTransitionProbabilityMats] = GET_TRANS_PROBS(partition(PL_start:PL_stop),subjInd(PL_start:PL_stop));
[PLTransitionProbabilityNoPersist2D,PLTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(partition(PL_start:PL_stop),subjInd(PL_start:PL_stop));


%Overall
grpAvgLSD = squeeze(mean(LSDTransitionProbabilityMats,1));% .* ~eye(numClusters);
% % nans occur for transitions from states that are not present at all for a subject or within the tested blocks. 
% % transitions to that state are 0
grpAvgPL = squeeze(mean(PLTransitionProbabilityMats,1));% .* ~eye(numClusters); 


%% permutation testing or t-test to compare transition probability matrices
% %average
[~,pavg,~,t]=ttest(LSDTransitionProbabilityMats,PLTransitionProbabilityMats);
% pavg = triu(reshape(pavg,numClusters,numClusters),1) + tril(reshape(pavg,numClusters,numClusters),-1);
fdravg = mafdr(reshape(pavg,1,numClusters^2),'BHFDR',1);
fdravg = reshape(fdravg,numClusters,numClusters);

% fdravg(fdravg==0) = NaN; %do not run this line and the one that follows if you are plotting persist method
% pavg(pavg==0) = NaN; %will use this for plotting non-corrected values




%% plot transition probabilities
maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales

f = figure;

subplot(1,3,1);
imagesc(grpAvgLSD);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvgPL);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);

ylabel('Current State'); xlabel('Next New State');
title('PL');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,3);
LSDMinusPLTP = squeeze(t.tstat);%(grpAvgPL - grpAvgLSD); %switching order for manuscript figures!
% imagesc(LSDMinusPLTP.*~eye(numClusters)); colormap('plasma'); %choose between this line and next line based on whether you are running NoPersist or regular TP 
imagesc(LSDMinusPLTP); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Next State');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y,x] = find(fdravg < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x-.12,y+.12,'**','Color','w','Fontsize', 32);
[y2,x2] = find(squeeze(pavg) < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x2-.12,y2+.12,'*','Color','w','Fontsize', 32);
caxis_bound = 4;%max(max(abs(LSDMinusPLTP)));
h = colorbar; ylabel(h,'t-stat'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('PL - LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');


%check [x,y] inputs ins plot 3 for corrected vs non-corrected pvalue


%%
clusters=char(clusterNames);
save(fullfile(savedir,['TransProbsData_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'LSDTransitionProbabilityMats','PLTransitionProbabilityMats','clusters','LSDTransitionProbabilityMatsNoPersist','PLTransitionProbabilityMatsNoPersist','LSDTransitionProbability2D','LSDTransitionProbabilityNoPersist2D','PLTransitionProbability2D','PLTransitionProbabilityNoPersist2D');
