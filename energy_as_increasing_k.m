%% SI Figure 14c make a plot to show the general claim that LSD flattens the landscape regarless of parameterization 

%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
%% set inputs

split='main'
load(fullfile(['data/',split,'.mat']))

distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000; % how many times to let kmeans try to converge before aborting rep


tot=nscans*2;

%% E vars


c = 0; T = 0.001 % set time scale parameters based on values from T_sweep_sps.m
Anorm = NORMALIZE(sc,c); 



%% load partitions from elbow plot

maxk = round(sqrt(TR)) -1; %k^2 must be less than TR to capture all transitions
k_rng=2:14;

C_sub = [repmat([1 0 0],nsubjs,1); repmat([0 0 0],nsubjs,1)];
C_avg = [[1 0 0];[0 0 0]];

figure()
set(gca,'FontSize',12);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','arial');
hold on
xlim([1 15]);
xlabel('Number of Clusters (k)');
ylabel('E_a_v_g');
for numClusters = k_rng
    load(fullfile(basedir,['/repkmeans/kmeans',num2str(split),'k_',num2str(numClusters),'.mat']));
    
    Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
    Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = [];
    
    LSDcentroids=NaN(nsubjs,nparc,numClusters);
    PLcentroids=NaN(nsubjs,nparc,numClusters);
    
    
    E_full_LSD=NaN(nsubjs,numClusters^2);
    E_full_PL=NaN(nsubjs,numClusters^2);
    
    for i=1:nsubjs
        LSDcentroids(i,:,:) = GET_CENTROIDS(concTS(LSDsubjInd==i,:),partition(LSDsubjInd==i),numClusters);
        PLcentroids(i,:,:) = GET_CENTROIDS(concTS(PLsubjInd==i,:),partition(PLsubjInd==i),numClusters);
        
         x0 = squeeze(LSDcentroids(i,:,Xo_ind));
         xf = squeeze(LSDcentroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
         WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
         E_full_LSD(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false);
         
         x0 = squeeze(PLcentroids(i,:,Xo_ind));
         xf = squeeze(PLcentroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
         WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
         E_full_PL(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false);
    end
    
    LSD = mean(E_full_LSD,2);
    PL = mean(E_full_PL,2);
    [~,pp(numClusters),~,t] = ttest(LSD,PL);
    
    scatter(repmat(numClusters,nsubjs*2,1),[LSD;PL],100,C_sub,'jitter','on','jitterAmount', 0.15,'LineWidth',1.5); 
    scatter([numClusters,numClusters],mean([LSD,PL]),200,C_avg,'filled');  
%     if pp(numClusters)<0.0001
%         text(numClusters,mean(mean([LSD,PL])),num2str(pp));
        text(numClusters,2950000,'*','Fontsize', 26);
%     end
end



%%
