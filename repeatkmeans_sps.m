clear all; close all;clc


%% set directories and load data

basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);

%load TxnParc matrix
load(fullfile(basedir,['LSD_ls463_cat.mat'])) % make sure this file matches the split you want to run below
addpath(genpath('code'))


savedir = fullfile(basedir,'repkmeans');
mkdir(savedir); cd(savedir);


%% set inputs

concTS = TS_ls463_nomean; %make sure this matches the split

split = 22; % i am using split to denote different processing applied to data 

%option to z-score time series
zdim = 0; 
nreps = 50; % how many times to repeat clustering. will choose lowest error solution
distanceMethod = 'correlation';
maxI = 1000; % how many times you allow kmeans to try to converge

if zdim > 0
	concTS = zscore(concTS,[],zdim);
end

disp('data loaded');

for numClusters = 2:14 %because we have 220 frames, k^2 must be less than 220 to capture all transitions
    
    disp(['K = ',num2str(numClusters),' Split = ',num2str(split)]);

%     savedir = ['kmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
%     mkdir(savedir);
%     cd(savedir);

    disp('start k-means');
    
    [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    save(fullfile(savedir,['kmeans',num2str(split),'k_',num2str(numClusters),'.mat']),'partition','sumd')
   
    clear partition
    clear sumd
    
    
end

disp('complete');
