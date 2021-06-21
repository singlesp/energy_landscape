%% generate subject specific centroids to make subj-specific comparisons and an TE comparison matrix

clear all; close all;
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

%% load BOLD data

% replace TS with your BOLD data formatted as a T-by-nparc matrix
% where T is the number of time points and nparc is the number of ROIs
load(fullfile(basedir,['LSD_ls463_cat.mat']),'TS_ls463_nomean') % make sure this file matches the split you want to run (see above)

split=22;

TR=217;
concTS = TS_ls463_nomean; %zscore(TS_ls463_nomean); %make sure this matches the split
[~,nparc] = size(concTS);
numClusters=4;
nsubjs=15;
nscans=29;
tot=nscans*2;
numNets=7;
D = NaN(nsubjs,TR*2,numClusters); %distance matrices will be stored here
LSDsubjInd=[repelem([1 3:nsubjs],TR),repelem(1:nsubjs,TR),repelem(0,TR*29)]'; % index LSD data from each subject
PLsubjInd=[repelem(0,TR*29),repelem([1 3:nsubjs],TR),repelem(1:nsubjs,TR)]';

% LSDsubjInd=[repelem(1:nsubjs,TR),repelem(0,TR*nsubjs)]'; % index data from each subject (how most data is organized)
% PLsubjInd=[repelem(0,TR*nsubjs),repelem(1:nsubjs,TR)]';

LSD_stop=TR*nscans;
PL_start=LSD_stop+1;
PL_stop=LSD_stop*2;
parts=NaN(TR*2,tot);

clusterNames={};
clusterNamesUp={};
clusterNamesDown={};


%%

%load partition to calculate centroids on

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'partition','centroids'); %make sure file matches centroids you want to initialize from

overallCentroids=centroids;

clear centroids;
centroids=NaN(nscans,nparc,numClusters);

for i=1:nsubjs
    
    disp(['Generating LSD centroids for subj ',num2str(i)]);
    centroids(i,:,:) = GET_CENTROIDS(concTS(LSDsubjInd==i,:),partition(LSDsubjInd==i),numClusters);
    [clusterNames{i},~,~,net7angle(i,:,:)] = NAME_CLUSTERS_ANGLE(squeeze(centroids(i,:,:)));
    [clusterNamesUp{i},clusterNamesDown{i},net7angle_Up(i,:,:),net7angle_Down(i,:,:)] = NAME_CLUSTERS_UP_DOWN(squeeze(centroids(i,:,:)));
   
    disp(['Generating PL centroids subj ',num2str(i)]);
    centroids(i+nsubjs,:,:) = GET_CENTROIDS(concTS(PLsubjInd==i,:),partition(PLsubjInd==i),numClusters);
    [clusterNames{i+nsubjs},~,~,net7angle(i+nsubjs,:,:)] = NAME_CLUSTERS_ANGLE(squeeze(centroids(i+nsubjs,:,:)));
    [clusterNamesUp{i+nsubjs},clusterNamesDown{i+nsubjs},net7angle_Up(i+nsubjs,:,:),net7angle_Down(i+nsubjs,:,:)] = NAME_CLUSTERS_UP_DOWN(squeeze(centroids(i+nsubjs,:,:)));
    
end

%% ttest on radial values


for i=1:numClusters
    [~,pUP(i,:),~,statUP{i}]=ttest(net7angle_Up(1:nsubjs,i,:),net7angle_Up(nsubjs+1:nsubjs*2,i,:));
    [~,pDOWN(i,:),~,statDOWN{i}]=ttest(net7angle_Down(1:nsubjs,i,:),net7angle_Down(nsubjs+1:nsubjs*2,i,:));
end

pfdrDN=reshape(mafdr(reshape(pDOWN,1,[]),'BHFDR',1),numClusters,numNets);
pfdrUP=reshape(mafdr(reshape(pUP,1,[]),'BHFDR',1),numClusters,numNets);


UPstat = NaN(numClusters,numNets);
DNstat = NaN(numClusters,numNets);
for i=1:numClusters
    UPstat(i,:) = squeeze(statUP{i}.tstat(:,:,1:7));
    DNstat(i,:) = squeeze(statDOWN{i}.tstat(:,:,1:7));
end

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
clusterNames = {'MS-1a','MS-1b','MS-2a','MS-2b'};

figure;
subplot(1,2,1);
imagesc(UPstat); colormap('plasma');
xticks(1:numNets); xticklabels(YeoNetNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Centroids'); xlabel('Networks');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y1,x1] = find(pUP < sig_thresh); 
text(x1-.12,y1+.12,'*','Color','w','Fontsize', 24);
[y,x] = find(pfdrUP < sig_thresh); 
text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
u_caxis_bound = max(max(4));
l_caxis_bound = min(min(-4));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('UP Radial values from subj centroids LSD > PL');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');


subplot(1,2,2);
imagesc(DNstat); colormap('plasma');
xticks(1:numNets); xticklabels(YeoNetNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Centroids'); xlabel('Networks');
sig_thresh = 0.05; %/numClusters^2; % bonferroni correction, for two-tailed p-values so only
[y1,x1] = find(pDOWN < sig_thresh); 
text(x1-.12,y1+.12,'*','Color','w','Fontsize', 24);
[y,x] = find(pfdrDN < sig_thresh); %this was used with two-tailed perm test-> p___.*~eye(numClusters)
text(x-.12,y+.12,'**','Color','w','Fontsize', 24);
u_caxis_bound = max(max(4));
l_caxis_bound = min(min(-4));
h = colorbar; ylabel(h,'t-stat'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('DOWN Radial values from subj centroids LSD > PL');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

% save sch454radialvals.mat net7angle_Up net7angle_Down

%% make a correlation plot comparing LSD average centroids with PL average centriods
load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'clusterNames'); 
LSDcentroids = squeeze(mean(centroids(1:nsubjs,:,:)));
PLcentroids = squeeze(mean(centroids(nsubjs+1:nsubjs*2,:,:)));

err = NaN(numClusters,numClusters);
for i=1:numClusters
    for k=1:numClusters
        err(i,k)=immse(LSDcentroids(:,i),PLcentroids(:,k));
    end
end

[rcent,pcent]=corr(LSDcentroids,PLcentroids);

figure;
imagesc(rcent);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('LSD Centroids'); xlabel('PL Centroids');
title('Correlation between Condition-Specific Centroids');
caxis_bound = 1;
h = colorbar; ylabel(h,'R'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

figure;
imagesc(err);
xticks(1:numClusters); yticks(1:numClusters); colormap('hot'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('LSD Centroids'); xlabel('PL Centroids');
title('MSE between Condition-Specific Centroids');
u_caxis_bound = min(max(err));
l_caxis_bound = min(min(err));
h = colorbar; ylabel(h,'MSE'); caxis([l_caxis_bound u_caxis_bound]); h.Ticks = [l_caxis_bound u_caxis_bound]; h.TickLabels = [round(l_caxis_bound,2,'significant') round(u_caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');


%% plot radial plots for each subj/condition (can skip unless interested)

% YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};
% YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
% numNets = numel(YeoNetNames);
% clusterColors = GET_CLUSTER_COLORS(numClusters);
% clusterColors = hex2rgb(clusterColors);
% %YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;247 253 205; 218 166 86; 199 109 117] / 255;
% YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117] / 255;
% YeoColors = [YeoColors;YeoColors];
% netAngle = linspace(0,2*pi,numNets+1);
% thetaNames = YeoNetNames; thetaNames{8} = '';
% 
% 
% for i=1:nsubjs*2
% 
%     overallNames = clusterNames{i};
%     icentroids = squeeze(centroids(i,:,:));
%     inet7angle = squeeze(net7angle(i,:,:));
%     inet7angle_Up = squeeze(net7angle_Up(i,:,:));
%     inet7angle_Down = squeeze(net7angle_Down(i,:,:));
%     f(i)=figure;
%     for K = 1:numClusters
%         ax = subplot(1,numClusters,K,polaraxes); hold on
%         polarplot(netAngle,[inet7angle_Up(K,:) inet7angle_Up(K,1)],'k');
%         polarplot(netAngle,[inet7angle_Down(K,:) inet7angle_Down(K,1)],'r');
%         thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
%         rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
%         for L = 1:numNets
%             ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%                 YeoColors(L,:), ax.ThetaTickLabel{L});
%         end
%         set(ax,'FontSize',10);
%         title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
%     end
%     f(i).PaperUnits = 'inches';
%     f(i).PaperSize = [8 1.5];
%     f(i).PaperPosition = [0 0 8 1.5];
%     
% end


%% energy matrices


if nparc == 454
    load Schaefer454_HCP_DTI_count.mat connectivity
%     load sch454_DTI_fiber_consensus_HCP.mat connectivity %consensus matrix
    
elseif nparc == 232
    load Schaefer232_HCP_DTI_count.mat connectivity

elseif nparc == 463
    load Lausanne463_HCP_DTI_count.mat connectivity
    
elseif nparc == 461
    load Lausanne463_HCP_DTI_count.mat connectivity
%     load ls463_DTI_fiber_consensus_HCP.mat connectivity %consensus matrix (less sparse)
    
    connectivity([14 463],:)=[];
    connectivity(:,[14 463])=[];
end

sc = connectivity;
c = 0; T = 0.001; % set time scale parameters based on values from T_sweep_sps.m
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
elseif nparc == 463
    load 5HTvecs_ls463.mat mean5HT2A_ls463
    HT = mean5HT2A_ls463;
elseif nparc == 461
    load 5HTvecs_ls463.mat mean5HT2A_ls463
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

%% plot either LSD v PL or PL weighted vs PL full

load(['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat'],'clusterNames'); 

% Energy = E_full;
Energy1=E_weighted(nsubjs+1:nsubjs*2,:);
Energy2=E_full(nsubjs+1:nsubjs*2,:);

% [~,pavg,~,t]=ttest(Energy(1:nsubjs,:),Energy(nsubjs+1:nsubjs*2,:));
[~,pavg,~,t]=ttest(Energy1,Energy2);
fdravg = mafdr(pavg,'BHFDR',1);
fdravg = reshape(fdravg,[numClusters numClusters])';
pavg = reshape(pavg,[numClusters numClusters])';

% grpAvgLSD = reshape(mean(Energy(1:nsubjs,:)),[numClusters numClusters])';
% grpAvgPL = reshape(mean(Energy(nsubjs+1:nsubjs*2,:)),[numClusters numClusters])';
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



%% compute transition energy for each state with average Placebo centroids


E_weighted_PLavg = zeros(1,numClusters^2);

x0 = squeeze(mean(centroids(16:30,:,Xo_ind)));
xf = squeeze(mean(centroids(16:30,:,Xf_ind))); 

for transition = 1:numClusters^2
    [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
    E_weighted_PLavg(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
end

%% save

% save('sch454_subjcentroids.mat','centroids','concTS','f','partition','overallCentroids','LSDsubjInd','PLsubjInd','MS1','MS1a','MS1b','MS2','MS2a','MS2b','E_full','E_weighted');
save(fullfile(savedir,['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat']),'centroids','E_weighted','E_full','T','E_weighted_PLavg')
