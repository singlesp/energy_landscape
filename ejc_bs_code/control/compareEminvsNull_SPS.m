%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 4;
split=22; % i am using split to denote different processing applied to data 
% split 0 corresponds to gsr only (LSDgsr_cat.mat)
% split 1 corresponds to gsr + bp 0.008-0.09 (LSDgsr_bp_cat.mat)
% split 2 corresponds to gsr + bp 0.008-0.2 (LSDgsr_bp2_cat.mat)
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
% nsubjs=15;
% nscans=29;
% TR=217;
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
c = 0; T = 0.001; % set time scale parameters based on values from paper
Anorm = NORMALIZE(sc,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

%% control energy -> full file here /Users/sps253/Documents/brain_states-master/code/control/transitionEnergyDynamicsGroupv2_TSweep_CFN.m

% centroids = squeeze(mean(centroids2(16:30,:,:)));% - 0.33 - randn(454,4);
% centroids(centroids<0) = 0;
% centroids(centroids>0) = 0;
% centroids = centroids +1;
% T = 0.1;
% centroids = zscore(centroids);

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

x0 = centroids(:,Xo_ind);
xf = centroids(:,Xf_ind); % now each column of x0 and xf represent state transitions
WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
E_full = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition

% compute weighted control energy:
% load('/Users/sps253/Documents/ROI_maps/sch454_to_yeo.csv'); 
% network7labels=sch454_to_yeo;
% InputVector = ismember(network7labels(1:nparc),2); % weight input towards a given system as in task

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
end

norm = (HT/max(HT))'; %zscore(HT); %
InputVector = norm; %(norm+3)/max(norm+3);% > 0.85; %option to binarize input vector
B = InputVector .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector
E_weighted = zeros(1,numClusters^2);
for transition = 1:numClusters^2    
    [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
    E_weighted(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
end

%% plot
% 
% LSD_E = mean(E_full(1:15,:));
% PL_E = mean(E_full(16:30,:));

Energy = (E_weighted_PLavg);
f=figure;
imagesc(reshape(Energy, [numClusters numClusters])'); colormap('viridis');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Initial State'); xlabel('Final State');
up_caxis_bound =  max(max(Energy));%.*~eye(numClusters));%max(max(abs(Energy)));
l_caxis_bound = min(min(Energy(offDiag)));%min(min(abs(Energy)));
h = colorbar; ylabel(h,'Emin (a.u.)'); caxis([l_caxis_bound up_caxis_bound]); %h.Ticks = [l_caxis_bound up_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(up_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('5-HT_2_a  > Uniform Inputs');% - Reshuffled E\_min');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

%% box plots

figure; 
hold on
boxplot(reshape(E_full, [numClusters numClusters])','colors','r');
boxplot(reshape(E_weighted, [numClusters numClusters])','colors','b');
ylim([l_caxis_bound-10000 up_caxis_bound+10000]);
title('Overall Energies (Columns)'); ylabel('Energy (a.u.)'); xticklabels(clusterNames);
% text(3,800,'Red = Uniformly-Weighted (a)','color','r');
% text(3,250,'Blue = 5-HT_2_a-Weighted (b)','color','b');
set(gca,'FontSize',18); set(gca,'TickLength',[0 0]); set(gca,'Fontname','arial');

%% compute transition energy for each state

NullTransitionEnergy = zeros(nperms,numClusters^2);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
% 	load(fullfile(datadir,'Group_DLWNull_Gramians',['DLWNull_FA_Laus',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'A_DLWnull');
	
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

E_matrix = reshape(E_weighted, [numClusters numClusters])';

% % reshape BCT null into matrix
% [num_transitions,nperms,~] = size(E_BCTnull);
% BCTNullTransitionEnergyMatrix = zeros(numClusters,numClusters,nperms);
% 
% for P = 1:nperms
% 	BCTNullTransitionEnergyMatrix(:,:,P) = reshape(E_BCTnull(:,P,T_idx),[numClusters numClusters])';
% end

% reshape DLW null into matrix
[nperms,num_transitions] = size(NullTransitionEnergy);
NullTransitionEnergyMatrix = zeros(numClusters,numClusters,nperms);

for P = 1:nperms
    NullTransitionEnergyMatrix(:,:,P) = reshape(NullTransitionEnergy(P,:),[numClusters numClusters])';
end

sig_thresh = 0.05 / num_transitions; % bonferroni correct over all transitions
pvals_onetail = mean(E_matrix > NullTransitionEnergyMatrix,3); % get p-value as % of times real energy is higher than energy in null matrix for a transition
% pvals_onetail_BCT = mean(E_matrix > BCTNullTransitionEnergyMatrix,3); % get p-value as % of times real energy is higher than energy in null matrix for a transition

%%
Energydiff = reshape((E_weighted - mean(NullTransitionEnergy,1)),4,4)';


    
[h(1),p(1),~,~]=ttest(Energydiff(:,1),Energydiff(:,4));
[h(2),p(2),~,~]=ttest(Energydiff(:,1),Energydiff(:,3));
[h(3),p(3),~,~]=ttest(Energydiff(:,1),Energydiff(:,2));

[h(4),p(4),~,~]=ttest(Energydiff(:,2),Energydiff(:,4));
[h(5),p(5),~,~]=ttest(Energydiff(:,2),Energydiff(:,3));

[h(6),p(6),~,~]=ttest(Energydiff(:,3),Energydiff(:,4));

pfdr=mafdr(p,'BH',true);

%% plot energy matrix
f=figure;
imagesc(E_matrix); %reshape(Energy,[numClusters numClusters])');%E_matrix); 
h = colorbar; 
ylabel(h,'E_{min}');
xticks(1:numClusters); yticks(1:numClusters); colormap('viridis');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
[y,x] = find(pvals_onetail < sig_thresh);
text(x-.15,y+.15,'**','Color','w','Fontsize', 36);
[y2,x2] = find(pvals_onetail < 0.05); 
text(x2-.15,y2+.15,'*','Color','w','Fontsize', 36);
h = colorbar; ylabel(h,'Emin (a.u.)'); caxis([l_caxis_bound up_caxis_bound]);
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('5-HT_2_a Receptor-Weighted Inputs');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 6 6];
f.PaperSize = [6 6];

%% save
% 
clusters=char(clusterNames);
save(fullfile(savedir,['EnergyData_bp',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'full2Avec.mat']), 'E_full','E_weighted','clusters','NullTransitionEnergy','NullTransitionEnergyMatrix');
