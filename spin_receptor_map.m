%% spin test (spin receptor map - maps obtained from Vasa code) - Figure 4a/b iv
% ref: Vasa et al Cerebral Cortex 2018
% (https://doi.org/10.1093/cercor/bhx249)

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%%
split='main' 
load(fullfile(['data/',split,'.mat']))

numClusters=4; T=0.001; 


load(['subjenergy_split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat'], 'centroids', 'E_weighted_PLavg');
load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'], 'clusterNames');





%% load SC - for spin test to work, need to add missing region back in

if nparc == 461
    load data/Lausanne463_HCP_DTI_count.mat vol_normalized_sc% connectivity
%     load ls463_DTI_fiber_consensus_HCP.mat connectivity %consensus matrix (less sparse)

    A = centroids(:,1:13,:);
    B = centroids(:,244,:);
    C = centroids(:,14:end,:);
%     D = zeros(30,1,4); sans brainstem
    
    centroids = horzcat(A,B,C);
    
    nparc=462;
    
    connectivity = vol_normalized_sc;
    connectivity(463,:)=[]; %sans brainstem
    connectivity(:,463)=[];
    sc = connectivity;
end

%%


c = 0; 
Anorm = NORMALIZE(sc,c); 

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities
if nparc == 454
    load data/5HTvecs_sch454.mat mean5HT2A_sch454
    HT = mean5HT2A_sch454;
elseif nparc == 232
    load data/5HTvecs_sch232.mat mean5HT2A_sch232
    HT = mean5HT2A_sch232;

elseif nparc == 462
    load data/5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
    

end %weight towards 5HT2a
norm_HT = (HT/max(HT))';

%% spin test

if nparc ==232
    load(['data/SpinTests/rotated_maps/rotated_Schaefer_200.mat'])
    
    load('data/sch232_to_yeo.csv');
    network7labels=sch232_to_yeo';
    nosub_HT=norm_HT;
    nosub_HT(network7labels>7)=[];

    
    
    spun_2a = [];
    
    for i=1:size(perm_id,2)
       
        for r=1:size(perm_id,1)
            
            spun_2a(r,i) = nosub_HT(perm_id(r,i));
            
        end
    end
        
    a = spun_2a; % cortical assignments
    b = reshape(repelem(norm_HT(201:232),10000),10000,32)'; % subcortical assignments
  
    perm_map_2a = vertcat(a,b);
    
      
    
end

if nparc ==462
    load(['data/SpinTests/rotated_maps/rotated_Lausanne_463.mat'])
    
    load('data/Lausanne_463_subnetworks.mat');
    network7labels=subnetworks_reorder;
    nosub_HT=norm_HT;
    nosub_HT(network7labels>7)=[];

    perm_id = perm_id+1;
    
    spun_2a = [];
    
    for i=1:size(perm_id,2)
       
        for r=1:size(perm_id,1)
            
            spun_2a(r,i) = nosub_HT(perm_id(r,i));
            
        end
    end
        
    a = spun_2a(1:223,:); %right hemi cortical assignments
    b = reshape(repelem(norm_HT(224:230),10000),10000,7)'; %right hemi subcortical assignments
    c = spun_2a(224:end,:); %left hemi cortical assignments
    d = reshape(repelem(norm_HT(456:463),10000),10000,8)'; %left hemi subcortical assignments
    
    perm_map_2a = vertcat(a,b,c,d);
    
    perm_map_2a(463,:) = []; %sans brainstem
    
end


%%
nperms=size(perm_map_2a,2);

x0 = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xo_ind)));
xf = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xf_ind))); 



NullTransitionEnergy = zeros(nperms,numClusters^2);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
	
    InputVector = perm_map_2a(:,P)'; %randomly reshuffle input vector
    B = InputVector .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector
    E_weighted_RAND = zeros(1,numClusters^2);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_weighted_RAND(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
	NullTransitionEnergy(P,:) = E_weighted_RAND;
end

%% compute pvalues

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

clusterNames={'MS-1a';'MS-1b';'MS-2a';'MS-2b'};

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
u_caxis_bound = .5;%max(max(LSDMinusPLTP));
l_caxis_bound = 0;%min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'corrected p'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,1,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('Shuffled > True');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

%% save
save(fullfile(savedir,['receptorSPIN__split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat']),'E_weighted_PLavg','NullTransitionEnergy','pvals_onetail','fdrpv1t','nperms');