%% SI Figure 2

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%%

split='main' 
load(fullfile(['data/',split,'.mat']))
numClusters = 4;

nobs = nsubjs


%% create null and random TS

for N = 1:nobs
	disp(['Subject ', num2str(N)]);
    % rest independent phase randomization
    concTS(LSDsubjInd == N,:) = linsurr_ind(concTS(LSDsubjInd == N,:));
    % nback independent phase randomization
    concTS(PLsubjInd == N,:) = linsurr_ind(concTS(PLsubjInd == N,:));    
end
iprTS = concTS;

% save(fullfile(savedir,['iprtimeseries',num2str(split),'.mat']),'iprTS');
clear concTS

randTS = randn(size(iprTS));
% save(fullfile(savedir,['randtimeseries',num2str(split),'.mat']),'randTS');

%% cluster null and random TS
numClusters = 4;
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000;

[partition_ipr] = kmeans(iprTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
% save(fullfile(savedir,['IPRpartition_k',num2str(numClusters)]),'partition_ipr')

[partition_randn] = kmeans(randTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);;
% save(fullfile(savedir,['randnpartition_k',num2str(numClusters)]),'partition_randn')

%%


n_sil = 10;
% sil_mask = ismember(subjInd,randperm(nobs,n_sil));

load(fullfile(['data/',split,'.mat']),'concTS')

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'partition');

f = figure; 
subplot(1,4,1);
[S_randn,~] = silhouette(randTS,partition_randn,distanceMethod);
title({'Independent','Random Gaussians'});
set(gca,'FontSize',8);

subplot(1,4,2);
[S_ipr,~] = silhouette(iprTS,partition_ipr,distanceMethod);
title({'Autocorrelation-Preserving','Null'});
set(gca,'FontSize',8);

subplot(1,4,3);
[S_data,~] = silhouette(concTS,partition,distanceMethod);
title('Real Data');
set(gca,'FontSize',8);

subplot(1,4,4);
histogram(S_data,'EdgeAlpha',0.3); hold on;
histogram(S_ipr,'EdgeAlpha',0.3);
xlim([min([S_data;S_ipr]) max([S_data;S_ipr])]);
xlabel('Silhouette Value');
ylabel('# of TRs');
[~,p,ci] = ttest2(S_data,S_ipr);
title(['\mu_{real}-\mu_{ipr} = ',num2str(round(mean(ci),2,'significant')),', p = ',num2str(round(p,2,'significant'))])
set(gca,'FontSize',8);

% f.PaperUnits = 'centimeters';
% f.PaperSize = [21 4];
% f.PaperPosition = [0 0 21 4];
% set(0,'CurrentFigure',f);
% saveas(f,fullfile(savedir,['SilhouetteVsNulls_k',num2str(numClusters),'.pdf']));
%print(fullfile(savedir,'SilhouetteVsNulls.png'),'-dpng','-r400');

%%

f = figure;
histogram(S_ipr); hold on;
histogram(S_data);
xlabel('Silhouette Value');
ylabel('Count');
[~,p,ci,stats] = ttest2(S_data,S_ipr);
stats.tstat
stats.df
title(['\mu_{real}-\mu_{ipr} = ',num2str(round(mean(ci),2,'significant')),', p = ',num2str(round(p,2,'significant'))])
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 8 8];
f.PaperSize = [8 8];
% saveas(f,fullfile(savedir,['SilhouetteRealvsIPRHistogram.pdf']));