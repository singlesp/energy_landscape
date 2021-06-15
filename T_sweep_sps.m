%% find optimum T for LSD and PL
clear all; close all;
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);



%% set inputs
numClusters = 4;
split=22; % i am using split to denote different processing applied to data 

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

load(fullfile(savedir,['TransProbsData_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
%%%%% LOAD SUBJECT SPECIFIC CENTROIDS HERE. MAKE SURE YOU'RE LOADING FROM
%%%%% THE RIGHT PARC/SPLIT
load(fullfile(savedir,['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat']),'centroids')
LSD_centroids = squeeze(mean(centroids(1:15,:,:)));
PL_centroids = squeeze(mean(centroids(16:30,:,:)));

load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'centroids')
[nparc,~] = size(centroids);

nperms=1000;

% load sc_fc.mat sc90_nonsymm % load example group average structural A matrix -- this is from PNC while fMRI is HCP so DTI is younger than fMRI here
if nparc == 454
    load Schaefer454_HCP_DTI_count.mat connectivity
%     load sch454_DTI_fiber_consensus_HCP.mat connectivity %consensus matrix
    
elseif nparc == 232
    load Schaefer232_HCP_DTI_count.mat connectivity
    
elseif nparc == 461
    load Lausanne463_HCP_DTI_count.mat connectivity
%     load ls463_DTI_fiber_consensus_HCP.mat connectivity %consensus matrix (less sparse)
    
    connectivity([14 463],:)=[];
    connectivity(:,[14 463])=[];
end

sc = connectivity;

c = 0; % set time scale parameters based on values from paper
Anorm = NORMALIZE(sc,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

%% Compare group average E_full with Transprobs

T_rng = [0.001:0.5:10]; nT = length(T_rng);

Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

x0 = centroids(:,Xo_ind);
xf = centroids(:,Xf_ind);

x0_L = LSD_centroids(:,Xo_ind);
xf_L = LSD_centroids(:,Xf_ind);

x0_P = PL_centroids(:,Xo_ind);
xf_P = PL_centroids(:,Xf_ind);
% now each column of x0 and xf represent state transitions

% load 5HTvecs_sch454.mat mean5HT2A_sch454
% HT = mean5HT2A_sch454;
% norm = (HT/max(HT))'; %zscore(HT); %
% InputVector = norm; %(norm+3)/max(norm+3);% > 0.85; %option to binarize input vector
% B = InputVector .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector



  
E_full_grpavg_T = NaN(nT,numClusters^2);
E_full_LSD_T = NaN(nT,numClusters^2);
E_full_PL_T = NaN(nT,numClusters^2);

E_w_grpavg_T = NaN(nT,numClusters^2);
E_w_LSD_T = NaN(nT,numClusters^2);
E_w_PL_T = NaN(nT,numClusters^2);
  
for i=1:nT
    T=T_rng(i);

    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full_grpavg_T(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    E_full_LSD_T(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_L, xf_L, T,false);
    E_full_PL_T(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_P, xf_P, T,false);
    
%     for transition = 1:numClusters^2    
%         [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
%         E_weighted(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
%     end

end

%% corrs Persist

for i=1:nT
    
    [Rgrp_LTP(i),Pgrp_LTP(i)] = corr(E_full_grpavg_T(i,:)',mean(LSDTransitionProbability2D)','type','Spearman');
    [Rgrp_PTP(i),Pgrp_PTP(i)] = corr(E_full_grpavg_T(i,:)',mean(PLTransitionProbability2D)','type','Spearman');
    
    [R_L(i),P_L(i)] = corr(E_full_LSD_T(i,:)',mean(LSDTransitionProbability2D)','type','Spearman');
    
    [R_P(i),P_P(i)] = corr(E_full_PL_T(i,:)',mean(PLTransitionProbability2D)','type','Spearman');
end

%% corrs NoPersist

% for i=1:nT
%     
%     [Rgrp_LTP(i),Pgrp_LTP(i)] = corr(E_full_grpavg_T(i,offDiag)',mean(LSDTransitionProbability2D(:,offDiag))');%,'type','Spearman');
%     [Rgrp_PTP(i),Pgrp_PTP(i)] = corr(E_full_grpavg_T(i,offDiag)',mean(PLTransitionProbability2D(:,offDiag))');%,'type','Spearman');
%     
%     [R_L(i),P_L(i)] = corr(E_full_LSD_T(i,offDiag)',mean(LSDTransitionProbability2D(:,offDiag))');%,'type','Spearman');
%     
%     [R_P(i),P_P(i)] = corr(E_full_PL_T(i,offDiag)',mean(PLTransitionProbability2D(:,offDiag))');%,'type','Spearman');
% end

%% plot

rank2=tiedrank(mean(LSDTransitionProbability2D));
rank3=tiedrank(mean(PLTransitionProbability2D));

for i=1:11
    rank1=tiedrank(E_full_grpavg_T(i,:));
    rank4=tiedrank(E_full_LSD_T(i,:));
    rank5=tiedrank(E_full_PL_T(i,:));
    
    figure;
    subplot(2,2,1)
    scatter(rank1,rank2); lsline; title(['E\_full grp avg T=',num2str(T_rng(i)),' vs LSD Transition Probability']);
    xlabel('Energy (rank)'); ylabel('Transition Probability (rank)');
    text(12,14,['r = ',num2str(Rgrp_LTP(i))]); text(12,13,['p = ',num2str(Pgrp_LTP(i))]);
        
    subplot(2,2,2)
    scatter(rank1,rank3); lsline; title(['E\_full grp avg T=',num2str(T_rng(i)),' vs PL Transition Probability']);
    xlabel('Energy (rank)'); ylabel('Transition Probability (rank)');
    text(12,14,['r = ',num2str(Rgrp_PTP(i))]); text(12,13,['p = ',num2str(Pgrp_PTP(i))]);
    
    subplot(2,2,3)
    scatter(rank4,rank2); lsline; title(['E\_full LSD avg T=',num2str(T_rng(i)),' vs LSD Transition Probability']);
    xlabel('Energy (rank)'); ylabel('Transition Probability (rank)');
    text(12,14,['r = ',num2str(R_L(i))]); text(12,13,['p = ',num2str(P_L(i))]);
    
    subplot(2,2,4)
    scatter(rank5,rank3); lsline; title(['E\_full PL avg T=',num2str(T_rng(i)),' vs PL Transition Probability']);
    xlabel('Energy (rank)'); ylabel('Transition Probability (rank)');
    text(12,14,['r = ',num2str(R_P(i))]); text(12,13,['p = ',num2str(P_P(i))]);

end

%%
save(fullfile(savedir,['T_sweep_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'E_full_grpavg_T','E_full_LSD_T','E_full_PL_T');
