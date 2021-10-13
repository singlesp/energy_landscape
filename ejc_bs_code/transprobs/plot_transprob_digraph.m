% makes Fig 4e-f
%%
% plot transition probability as directed graph
addpath(genpath(fullfile(basedir,'code')));

%% load trans probs

% R = load(['results/',name_root,'/analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',...
%     num2str(numClusters),name_root,'.mat']);
LSDTransProbsFull = reshape(mean(LSDTransitionProbabilityMatsNoPersist,1),[numClusters numClusters])';

% N = load(['results/',name_root,'/analyses/nbackblocks/TransProbsNoPersist2back_k',...
%     num2str(numClusters),name_root,'.mat']);
PLTransProbsFull = reshape(nanmean(PLTransitionProbabilityMatsNoPersist,1),[numClusters numClusters])';

%% plot

thresh = 0.25;
% reformat trans prob matrix so that each edge is unidirectional
LSDpprestTransProbs = round(LSDTransProbsFull,2,'significant');% - restTransProbsFull';
LSDpprestTransProbs(LSDpprestTransProbs < thresh) = 0;
R = digraph(LSDpprestTransProbs);

PLppTransProbs = round(PLTransProbsFull,2,'significant');% - nBackTransProbsFull';
PLppTransProbs(PLppTransProbs < thresh) = 0;
N = digraph(PLppTransProbs);

all_edges = [R.Edges.Weight;N.Edges.Weight];
max_edge = max(all_edges(all_edges >0));
min_edge = min(all_edges(all_edges >0));
RWidths = 0.1+5*((R.Edges.Weight-min_edge)/max_edge);
NWidths = 0.1+5*((N.Edges.Weight-min_edge)/max_edge);

f = figure; 
subplot(1,2,1); plot(R,'Layout','circle','LineWidth',RWidths,'EdgeColor',[0 .361 .6335]); axis off; % for text add 'EdgeLabel',R.Edges.Weight,
title('LSD');
set(gca,'FontSize',24);
subplot(1,2,2); plot(N,'Layout','circle','LineWidth',NWidths,'EdgeColor',[0 .361 .6335]); axis off % for text add 'EdgeLabel',N.Edges.Weight,
title('PL');
set(gca,'FontSize',24);
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 12 4];
f.PaperSize = [12 4];
%saveas(f,fullfile(basedir,'results',name_root,'analyses/transitionprobabilities',['TransProbDigraph_Thresh',num2str(thresh),'_k',num2str(numClusters),'.pdf']));

%% plot scale bars for width of lines
WidthTicks = [0.25 0.3 0.35]; % set gradations in width axis for scale bar
WidthScale = 0.1+5*((WidthTicks-min_edge)/max_edge);
f=figure; hold on;
for j = 1:length(WidthScale)
    line([1 2],[j j],'LineWidth',WidthScale(j),'Color','k')
    text(2.2,j,num2str(WidthTicks(j)),'FontSize',12);
end
xlim([0 3]); ylim([0 5]); axis off

f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 2 2];
f.PaperSize = [2 2];
%saveas(f,fullfile(basedir,'results',name_root,'analyses/transitionprobabilities',['TransProbDigraphWidthAxisScale_Thresh',num2str(thresh),'_k',num2str(numClusters),'.pdf']));
