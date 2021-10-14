%% SI Figure 15 iv: cluster individuals separately and compare AMI to original AMI and E

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%%
numClusters = 4;
split='main'
load(fullfile(['data/',split,'.mat']))


load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'partition','centroids')


distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000;



%% run k-means on each sub & get centroids

parts = {};

for i=1:nsubjs
    disp(['Starting sub: ',num2str(i)]);
    
    subTS = concTS(subjInd ==i,:);
    
    %two options: initialize on main analysis centroids (line below) - this
    %will make the centroid order meaningful OR random initialization. 
    
%     [parts{i},~,~,~] = kmeans(subTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxI,'Start',repmat(centroids',[1,1,nreps]));
    [parts{i},~,~,~] = kmeans(subTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxI);
   
    if i==2 && strcmp(split,'main') %if doing gsr or sch, would need to change that here
        LSDind = [repelem(1,TR),repelem(0,TR)];
    elseif strcmp(split,'main') %and here
        LSDind = [repelem(1,TR*2),repelem(0,TR*2)];
    else
        LSDind = [repelem(1,TR),repelem(0,TR)];
    end
    
    LSDcentroids(i,:,:) = GET_CENTROIDS(subTS(LSDind==1,:),parts{i}(LSDind==1),numClusters);
    PLcentroids(i,:,:) = GET_CENTROIDS(subTS(LSDind==0,:),parts{i}(LSDind==0),numClusters);
end

%% compare LSD vs PL E with individuals are clustered on their own

c = 0; T = 0.001%2 % set time scale parameters based on values from T_sweep_sps.m
Anorm = NORMALIZE(sc,c); 

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities


E_full_LSD=NaN(nsubjs,numClusters^2);
E_full_PL=NaN(nsubjs,numClusters^2);

for i=1:nsubjs
    x0 = squeeze(LSDcentroids(i,:,Xo_ind));
    xf = squeeze(LSDcentroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full_LSD(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
    
    x0 = squeeze(PLcentroids(i,:,Xo_ind));
    xf = squeeze(PLcentroids(i,:,Xf_ind)); % now each column of x0 and xf represent state transitions
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full_PL(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition

end

%% compare

clusterNames = {'MS-1a'; 'MS-1b'; 'MS-2a'; 'MS-2b'};

Energy1=E_full_LSD;
Energy2=E_full_PL;

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
title('PCB');
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
u_caxis_bound = -1 %max(max(LSDMinusPLTP));
l_caxis_bound = -8 %min(min(LSDMinusPLTP));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound]; 
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('LSD > PCB');
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

%% compare 

LSD = mean(E_full_LSD,2)
PL = mean(E_full_PL,2)
[~,pp,~,t] = ttest(LSD,PL)

figure; violin([LSD,PL]); ylabel('E_a_v_g'); 
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
