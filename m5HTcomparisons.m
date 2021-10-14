%% Figure 5a/b
clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);

%% set inputs
numClusters = 4;
split='psilo' 
load(fullfile(['data/',split,'.mat']))

T = 0.001  % set time scale parameters based on values from paper
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'clusterNames');


c = 0; 
Anorm = NORMALIZE(sc,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

%% load subject-specific centroids and energies from subcentroids.m

load(fullfile(savedir,['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat']));
load(['subjenergy_split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat'], 'centroids','E_full', 'E_weighted');

E_l = E_full(1:nsubjs,:); %LSD centroid E calcs un-weighted
E_pl = E_full(nsubjs+1:nsubjs*2,:); %placebo centroid E calcs un-weighted

% Use this if not thresholding
E_w2a = E_weighted(nsubjs+1:nsubjs*2,:); %%placebo centroid E calcs 5-HT2a weighted

centroids = centroids(nsubjs+1:nsubjs*2,:,:); %placebo centroids

E_w1a = NaN(nsubjs,numClusters^2); %store 1a E's here
E_w1b = NaN(nsubjs,numClusters^2); %etc
E_w4 = NaN(nsubjs,numClusters^2);
E_wT = NaN(nsubjs,numClusters^2);

norm = NaN(nparc,5); %will want to scale receptor distributions to 0->1

%% control energy -> full file here /Users/sps253/Documents/brain_states-master/code/control/transitionEnergyDynamicsGroupv2_TSweep_CFN.m


% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities


% compute weighted control energy:

if nparc == 454
    load 5HTvecs_sch454.mat mean5HT2A_sch454
    HT = mean5HT2A_sch454;
elseif nparc == 232
    load 5HTvecs_sch232.mat mean5HT2A_sch232
    HT = mean5HT2A_sch232;
elseif nparc == 462
    load 5HTvecs_ls463.mat
    HT = horzcat(mean5HT2A_ls463,mean5HT1A_ls463,mean5HT1B_ls463,mean5HT4_ls463,mean5HTT_ls463);
    HT(463,:)=[];
    for i =1:5
        norm(:,i) = (HT(:,i)/max(HT(:,i)));
    end
    
    % below will get you distributions that are level between the
        % different receptors - these results also hold
%     for i =1:5
%         norm(:,i) = tiedrank(-HT(:,i))./nparc;
%     end
elseif nparc == 461
    load 5HTvecs_ls463.mat
    HT = horzcat(mean5HT2A_ls463,mean5HT1A_ls463,mean5HT1B_ls463,mean5HT4_ls463,mean5HTT_ls463);
    HT([14 463],:)=[];
    
    %scale receptor distributions to have max value of 1
    for i =1:5
        norm(:,i) = (HT(:,i)/max(HT(:,i)));
    end
    
    % below will get you distributions that are level between the
        % different receptors - these results also hold
    %     for i =1:5
%         norm(:,i) = tiedrank(-HT(:,i))./nparc;
%     end

end

%% option: threshold to a common amount.

%Lasuanne has 448 cortical ROI's, making an average of 64 ROI/cortical
%network. In order to simulate the process of weighting towards a RSN as
%was done by Cornblath et. al. 2020, we will threshold each of our receptor
%maps to give a total of 64 regions for weighting

% thresh_range = [0.1:0.001:1]; n=length(thresh_range);
% thresh_norm = zeros(nparc,5);
% thresh = NaN(1,5);
% for k=1:n
%     for i=1:5
%         if sum(thresh_norm(:,i))==64
%             
%             if isnan(thresh(1,i))==1
%                 thresh(1,i) = thresh_range(k);
%             end
%             i=i+1;
%         else
%             thresh_norm(:,i) = norm(:,i)>thresh_range(k);
%         end
%     end
% end
%% calc energies

%toggle norm or thresh_norm - if calc thresh norm need to recalc E_w2a
% norm = thresh_norm;
% E_w2a = NaN(15,numClusters^2);

for i=1:nsubjs
    
    x0 = squeeze(centroids(i,:,Xo_ind));
    xf = squeeze(centroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    
    % compute weighted control energy:
    
    %5-HT2a
%     InputVector = norm(:,1)';
%     B = InputVector .*eye(nparc) + eye(nparc);
%     for transition = 1:numClusters^2
%         [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
%         E_w2a(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
%     end
    
    %5-HT1a
    InputVector = norm(:,2)';
    B = InputVector .*eye(nparc) + eye(nparc);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_w1a(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
    
    %5-HT1b
    InputVector = norm(:,3)';
    B = InputVector .*eye(nparc) + eye(nparc);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_w1b(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
    
    %5-HT4
    InputVector = norm(:,4)';
    B = InputVector .*eye(nparc) + eye(nparc);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_w4(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
    
    %5-HTT
    InputVector = norm(:,5)';
    B = InputVector .*eye(nparc) + eye(nparc);
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_wT(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end
    

end


%% calc State means (columns)

%column means for each subject and each weighting
for i=1:nsubjs
    E_l_i(i,:)=mean(reshape(E_l(i,:),[4 4])');
    E_pl_i(i,:)=mean(reshape(E_pl(i,:),[4 4])');
    E_w2a_i(i,:)=mean(reshape(E_w2a(i,:),[4 4])');
    E_w1a_i(i,:) = mean(reshape(E_w1a(i,:),[4 4])');
    E_w1b_i(i,:) = mean(reshape(E_w1b(i,:),[4 4])');
    E_w4_i(i,:) = mean(reshape(E_w4(i,:),[4 4])');
    E_wT_i(i,:) = mean(reshape(E_wT(i,:),[4 4])');
end

%overall means for each subject and each weighting
E_l_m = mean(E_l_i,2);
E_pl_m = mean(E_pl_i,2);
E_w2a_m = mean(E_w2a_i,2);
E_w1a_m = mean(E_w1a_i,2);
E_w1b_m = mean(E_w1b_i,2);
E_w4_m = mean(E_w4_i,2);
E_wT_m = mean(E_wT_i,2);

%data for box plots
data = horzcat(E_w2a_m,E_w1a_m,E_w1b_m,E_w4_m,E_wT_m,E_pl_m);
% data = horzcat(E_l_m,E_w2a_m,E_w1a_m,E_w1b_m,E_w4_m,E_wT_m,E_pl_m);

%ttest for boxplots
[h(1),p(1)] = ttest(E_w2a_m,E_w1a_m);
[h(2),p(2)] = ttest(E_w2a_m,E_w1b_m);
[h(3),p(3)] = ttest(E_w2a_m,E_w4_m);
[h(4),p(4)] = ttest(E_w2a_m,E_wT_m);
[h(5),p(5),~,t] = ttest(E_w2a_m,E_pl_m);

% [h(6),p(6)] = ttest(E_w2a_m,E_l_m);
% 
% [h(7),p(7)] = ttest(E_l_m,E_w1a_m);
% [h(8),p(8)] = ttest(E_l_m,E_w1b_m);
% [h(9),p(9)] = ttest(E_l_m,E_w4_m);
% [h(10),p(10)] = ttest(E_l_m,E_wT_m);
% [h(11),p(11)] = ttest(E_l_m,E_pl_m);

pfdr = mafdr(p,'BH',1);

%% box plots

receptorNames={'5-HT2a','5-HT1a','5-HT1b','5-HT4','5-HTT','Uniform'};
colors = vertcat([1 0 0], [0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]);

figure; 
boxplot(data,'colors',colors);

% ylim([l_caxis_bound-10000 up_caxis_bound+10000]);
title('Overall Energies (Columns)'); ylabel('Energy (a.u.)'); ylim([80000 2000000]);
xticklabels(receptorNames);  
% text(3,800,'Red = Uniformly-Weighted (a)','color','r');
% text(3,250,'Blue = 5-HT_2_a-Weighted (b)','color','b');
set(gca,'FontSize',18); set(gca,'TickLength',[0 0]); set(gca,'Fontname','arial');

%%
%change name based on thresholding or not (use prefix <THRESHOLD_>
save(fullfile(savedir,['placebo_HT_compare_split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat']),'centroids','E_w2a','E_pl','E_l','T','E_w1a','E_w1b','E_w4','E_wT');