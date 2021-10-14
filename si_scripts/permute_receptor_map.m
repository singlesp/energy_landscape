%% perm test (randomly scramble receptor map)
% this was used in bioRxiv v1 and v2 - later chose to spin instead of
% shuffle

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%%
split='main'
load(fullfile(['data/',split,'.mat']))

numClusters=4;

load(['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat'], 'centroids');
load(['subjenergy_split',num2str(split),'_k',num2str(numClusters),'.mat'], 'E_weighted_PLavg');
load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'], 'clusterNames');





%%

c = 0; T = 0.001 % set time scale parameters based on values from T_sweep_sps.m
Anorm = NORMALIZE(sc,c); 

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities
if nparc == 454
    load 5HTvecs_sch454.mat mean5HT2A_sch454
    HT = mean5HT2A_sch454;
elseif nparc == 232
    load 5HTvecs_sch232.mat mean5HT2A_sch232
    HT = mean5HT2A_sch232;
elseif nparc == 461
    load 5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
    HT([14 463],:)=[];
elseif nparc == 462
    load 5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
    HT(463,:)=[];
end %weight towards 5HT2a
norm = (HT/max(HT))';
InputVector = norm;% > 0.883; %option to binarize input vector
B = InputVector .*eye(nparc) + eye(nparc);


%%
nperms=5 %10000;

x0 = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xo_ind)));
xf = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xf_ind))); 

%!!!if threshold, recalc an E_weighted_PLavg
% E_weighted_PLavg = zeros(1,numClusters^2);
% 
% for transition = 1:numClusters^2
%     [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
%     E_weighted_PLavg(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
% end

NullTransitionEnergy = zeros(nperms,numClusters^2);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
	
    IV = InputVector(randperm(length(InputVector))); %randomly reshuffle input vector
    B = IV .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector
    E_weighted_RAND = zeros(1,numClusters^2);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_weighted_RAND(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
	NullTransitionEnergy(P,:) = E_weighted_RAND;
end

%% compute pvalues
% [~,T_idx] = min(abs(T_rng-T)); % find index of T closest to desired T (will be the T_min identified at the end of transitionEnergyDynamicsGroupv2_TSweep.m)
% E_matrix = reshape(E_weighted(:,T_idx)',[numClusters numClusters])';

E_matrix = reshape(E_weighted_PLavg, [numClusters numClusters])';


[nperms,num_transitions] = size(NullTransitionEnergy);
NullTransitionEnergyMatrix = zeros(numClusters,numClusters,nperms);

for P = 1:nperms
    NullTransitionEnergyMatrix(:,:,P) = reshape(NullTransitionEnergy(P,:),[numClusters numClusters])';
end

% sig_thresh = 0.05 / num_transitions; % bonferroni correct over all transitions
pvals_onetail = mean(E_matrix > NullTransitionEnergyMatrix,3); % get p-value as % of times real energy is higher than energy in null matrix for a transition

fdrpv1t = mafdr(reshape(pvals_onetail,1,numClusters^2),'BHFDR',1);
fdrpv1t = reshape(fdrpv1t,numClusters,numClusters);

%% visualize

randdist = mean(NullTransitionEnergyMatrix,3);
truedist = reshape(E_weighted_PLavg, [numClusters numClusters])';

maxVal = max(max([randdist,truedist])); % sync color scales
minVal = min(min([randdist,truedist]));

figure;
subplot(1,3,1);
imagesc(randdist);
xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next State');
title('Random');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,2);
imagesc(truedist);
xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next State');
title('True');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,3);
LSDMinusPLTP = fdrpv1t;%(grpDiff); %((grpAvgPL - grpAvgLSD)); %switching order for manuscript figures
imagesc(LSDMinusPLTP); colormap('viridis');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Initial State'); xlabel('Final State');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y,x] = find(fdrpv1t < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
[y2,x2] = find(pvals_onetail < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x2-.15,y2+.18,'*','Color','w','Fontsize', 36);
u_caxis_bound = .1;%max(max(LSDMinusPLTP));
l_caxis_bound = 0;%min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'corrected p'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,1,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('Shuffled > True');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

%% save
save(fullfile(['receptorperm_split',num2str(split),'_k',num2str(numClusters),'.mat']),'E_weighted_PLavg','NullTransitionEnergy','pvals_onetail','fdrpv1t','nperms');