%  make elbow plot (Fig S2a-b) of within cluster / total variance explained by clusters at different k values
% assumes correlation distance

%% init

clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
savedir = fullfile(basedir,'results','example');mkdir(savedir);
load(fullfile(basedir,['LSD_ls463_cat.mat'])) % make sure this file matches the split you want to run below

%% inputs


concTS=TS_ls463_nomean(6294:end,:); %make sure this matches the split
N = size(concTS,1); % number of observations

split = 30; % i am using split to denote different processing applied to data 
% split 0 corresponds to gsr only (LSDgsr_cat.mat)
% split 1 corresponds to gsr + bp 0.008-0.09 (LSDgsr_bp_cat.mat)
% split 2 corresponds to gsr + bp 0.008-0.2 (LSDgsr_bp2_cat.mat)
% split 3 corresponds to LSD ONLY gsr + bp 0.008-0.2  (LSDgsr_bp2_LSD.mat)
% split 4 corresponds to PL ONLY gsr + bp 0.008-0.2  (LSDgsr_bp2_PL.mat)
% split 5 corresponds to the data fully preprocessed by ICL, parcellated
% into sch200 +32Tian, and GSR applied after parcellation
% (LSD_sch232_cat.mat TS_gsr_sch232)
% split 6 is same as 5 but only the LSD data (no PL)
% (TS_gsr_sch232(1:6293,:))
% split 7 is same as 5 but only the PL data (no LSD)
% (TS_gsr_sch232(6294:12586,:))
k_rng = 2:14;
VarianceExplained = zeros(length(k_rng),1);

for numClusters = k_rng
	disp(['K = ',num2str(numClusters)])
	load(fullfile(basedir,['/repkmeans/kmeans',num2str(split),'k_',num2str(numClusters),'.mat']));
	kClusterCentroids = GET_CENTROIDS(concTS,partition,numClusters);
	VarianceExplained(numClusters - 1) = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

save(fullfile(savedir,['VarianceExplained',num2str(split),'.mat']),'VarianceExplained','k_rng');

% Fig S2a-b
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['ElbowPlot',num2str(split),'.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(savedir,['GainInVarianceExplained',num2str(split),'.pdf']));