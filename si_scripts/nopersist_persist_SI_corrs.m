%% SI Figure 23: plot no persist imagesc and corr off diags

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);	

%% set inputs
numClusters = 4;


i=1; k=0;
split={'main'; 'gsr'; 'sch'; 'music'; 'psilo'};
clusterNames = {'MS-1a';'MS-1b';'MS-2a';'MS-2b'};
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];

f = figure;
for s=1:length(split)


    % nsubjs = 15%16
    
    
    % set save directory
    
    load(fullfile(savedir,['TransProbsData_bp',num2str(split{s}),'_k',num2str(numClusters),'.mat']))
    
    % isolate off diagonal from linearized transition probabilities
    
    
    grpAvgLSD = squeeze(mean(LSDTransitionProbabilityMats,1));% .* ~eye(numClusters);
    grpAvgPL = squeeze(mean(PLTransitionProbabilityMats,1));% .* ~eye(numClusters);
    
    [~,pavg,~,t]=ttest(LSDTransitionProbabilityMats,PLTransitionProbabilityMats);
    % pavg = triu(reshape(pavg,numClusters,numClusters),1) + tril(reshape(pavg,numClusters,numClusters),-1);
    fdravg = mafdr(reshape(pavg,1,numClusters^2),'BHFDR',1);
    fdravg = reshape(fdravg,numClusters,numClusters);
    
    % fdravg(fdravg==0) = NaN; %do not run this line and the one that follows if you are plotting persist method
    % pavg(pavg==0) = NaN; %will use this for plotting non-corrected values
    
    maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales
    
    
    
    
    subplot(5,3,1+k)
    LSDMinusPLTP = squeeze(t.tstat);%(grpAvgPL - grpAvgLSD); %switching order for manuscript figures!
    % imagesc(LSDMinusPLTP.*~eye(numClusters)); colormap('plasma');
    imagesc(LSDMinusPLTP); colormap('plasma');
    xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
    yticks(1:numClusters); yticklabels(clusterNames); axis square
    ylabel('Current State'); xlabel('Next State');
    sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
    [y,x] = find(fdravg < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
    text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
    [y2,x2] = find(squeeze(pavg) < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
    text(x2-.12,y2+.12,'*','Color','w','Fontsize', 24);
    caxis_bound = 4;%max(max(abs(LSDMinusPLTP)));
    h = colorbar; ylabel(h,'t-stat'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
    COLOR_TICK_LABELS(true,true,numClusters);
    %title([{'PL - LSD'};{'Persist'}]);
    set(gca,'FontSize',10);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    persistOD = LSDMinusPLTP(offDiag);
    
    grpAvgLSD = squeeze(mean(LSDTransitionProbabilityMatsNoPersist,1)) .* ~eye(numClusters);
    grpAvgPL = squeeze(mean(PLTransitionProbabilityMatsNoPersist,1)) .* ~eye(numClusters);
    
    [~,pavg,~,t]=ttest(LSDTransitionProbabilityMatsNoPersist,PLTransitionProbabilityMatsNoPersist);
    pavg = triu(reshape(pavg,numClusters,numClusters),1) + tril(reshape(pavg,numClusters,numClusters),-1);
    fdravg = mafdr(reshape(pavg,1,numClusters^2),'BHFDR',1);
    fdravg = reshape(fdravg,numClusters,numClusters);
    
    fdravg(fdravg==0) = NaN; %do not run this line and the one that follows if you are plotting persist method
    pavg(pavg==0) = NaN; %will use this for plotting non-corrected values
    
    maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales
    
    
    
    subplot(5,3,2+k)
    LSDMinusPLTP = squeeze(t.tstat);%(grpAvgPL - grpAvgLSD); %switching order for manuscript figures!
    % imagesc(LSDMinusPLTP.*~eye(numClusters)); colormap('plasma');
    imagesc(LSDMinusPLTP,'AlphaData',double(~isnan(LSDMinusPLTP))); colormap('plasma');
    xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
    yticks(1:numClusters); yticklabels(clusterNames); axis square
    ylabel('Current State'); xlabel('Next New State');
    sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
    [y,x] = find(fdravg < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
    text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
    [y2,x2] = find(squeeze(pavg) < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
    text(x2-.12,y2+.12,'*','Color','w','Fontsize', 24);
    caxis_bound = 4;%max(max(abs(LSDMinusPLTP)));
    h = colorbar; ylabel(h,'t-stat'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
    COLOR_TICK_LABELS(true,true,numClusters);
    %title([{'PL - LSD'};{'No Persist'}]);
    set(gca,'FontSize',10);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','arial');
    
    NOpersistOD = LSDMinusPLTP(offDiag);
    
    [r,p] = corr(persistOD',NOpersistOD','type','Spearman');
    rank1=tiedrank(persistOD);
    rank2=tiedrank(NOpersistOD);
    
    
    subplot(5,3,3+k)
    scatter(rank1,rank2); lsline; %title([{'Spearman Correlation Between Persist and'};{'No Persist TP Results (off diagonal)'}]);
    xlabel('Persist t-stat (rank)'); ylabel('No Persist t-stat (rank)');
    text(6,6,[{['r = ',num2str(r)]}; {['p = ',num2str(p)]}],'FontSize',10);
        
%     sgtitle(splitname{i});
    
    i=i+1;
    k=k+3;
    
end
