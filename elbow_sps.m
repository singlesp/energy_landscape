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


concTS=TS_ls463_nomean; %make sure this matches the split
N = size(concTS,1); % number of observations

split = 22; % i am using split to denote different processing applied to data 

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
