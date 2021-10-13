%% compute subject-specific energy matrices (Figure 4a/b ii-iii) after choosing T with T_sweep_sps.m
clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%% load BOLD data

c = 0;
T = 0.001 % set time scale parameters based on values from T_sweep_sps.m

numClusters=4;
split='main'
load(fullfile(['data/',split,'.mat']))

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'partition','centroids','clusterNames'); %make sure file matches centroids you want to initialize from

overallCentroids=centroids;

load(fullfile(savedir,['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat']),'centroids')



tot=nscans*2;

%% energy matrices



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
elseif nparc == 463
    load data/5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
elseif nparc == 461
    load data/5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
    HT([14 463],:)=[];
end %weight towards 5HT2a
norm = (HT/max(HT))';
InputVector = norm;% > 0.883; %option to binarize input vector
B = InputVector .*eye(nparc) + eye(nparc);

E_full=NaN(nsubjs*2,numClusters^2);
E_weighted=NaN(nsubjs*2,numClusters^2);

for i=1:nsubjs*2 
    x0 = squeeze(centroids(i,:,Xo_ind));
    xf = squeeze(centroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    
    % compute weighted control energy:
    
 % construct input matrix allowing input only into selected regions in InputVector
    
    for transition = 1:numClusters^2
        [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
        E_weighted(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
    end

end

%% plot PL weighted vs PL full


Energy1=E_weighted(nsubjs+1:nsubjs*2,:);
Energy2=E_full(nsubjs+1:nsubjs*2,:);


[~,pavg,~,t]=ttest(Energy1,Energy2);
fdravg = mafdr(pavg,'BHFDR',1);
fdravg = reshape(fdravg,[numClusters numClusters])';
pavg = reshape(pavg,[numClusters numClusters])';

grpAvgLSD = reshape(mean(Energy1),[numClusters numClusters])';
grpAvgPL = reshape(mean(Energy2),[numClusters numClusters])';

grpDiff = reshape(squeeze(t.tstat),[numClusters numClusters])';% .* -log(fdravg); (add sign() around tstat if want

maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales
minVal = min(min([grpAvgLSD,grpAvgPL]));

figure;
subplot(1,3,1);
imagesc(grpAvgLSD);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvgPL);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('PL');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,3);
LSDMinusPLTP = (grpDiff);%fdrpv1t; %((grpAvgPL - grpAvgLSD)); %switching order for manuscript figures
imagesc(LSDMinusPLTP); colormap('parula');%colormap('viridis');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Initial State'); xlabel('Final State');
sig_thresh = 0.05;
[y,x] = find(pavg < sig_thresh);
text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
[y,x] = find(fdravg < sig_thresh);
text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
u_caxis_bound = max(max(LSDMinusPLTP));
l_caxis_bound = min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('PL > LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

%% plot LSD vs PL 

Energy = E_full;


[~,pavg,~,t]=ttest(Energy(1:nsubjs,:),Energy(nsubjs+1:nsubjs*2,:));
fdravg = mafdr(pavg,'BHFDR',1);
fdravg = reshape(fdravg,[numClusters numClusters])';
pavg = reshape(pavg,[numClusters numClusters])';

grpAvgLSD = reshape(mean(Energy(1:nsubjs,:)),[numClusters numClusters])';
grpAvgPL = reshape(mean(Energy(nsubjs+1:nsubjs*2,:)),[numClusters numClusters])';


grpDiff = reshape(squeeze(t.tstat),[numClusters numClusters])';% .* -log(fdravg); (add sign() around tstat if want

maxVal = max(max([grpAvgLSD,grpAvgPL])); % sync color scales
minVal = min(min([grpAvgLSD,grpAvgPL]));

figure;
subplot(1,3,1);
imagesc(grpAvgLSD);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvgPL);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('PL');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([minVal maxVal]); colorbar

subplot(1,3,3);
LSDMinusPLTP = (grpDiff);%fdrpv1t; %((grpAvgPL - grpAvgLSD)); %switching order for manuscript figures
imagesc(LSDMinusPLTP); colormap('parula');%colormap('viridis');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Initial State'); xlabel('Final State');
sig_thresh = 0.05;
[y,x] = find(pavg < sig_thresh);
text(x-.15,y+.18,'*','Color','w','Fontsize', 36);
[y,x] = find(fdravg < sig_thresh);
text(x-.15,y+.18,'**','Color','w','Fontsize', 36);
u_caxis_bound = max(max(LSDMinusPLTP));
l_caxis_bound = min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('PL > LSD');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');



%% compute transition energy for each state with average Placebo centroids


E_weighted_PLavg = zeros(1,numClusters^2);

x0 = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xo_ind)));
xf = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,Xf_ind))); 

for transition = 1:numClusters^2
    [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
    E_weighted_PLavg(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
end

%%
save(fullfile(savedir,['subjenergy_split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat']),'E_weighted_PLavg','E_full','E_weighted','centroids')