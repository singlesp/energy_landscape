clear all; close all;clc
a=clock;
rng(a(6));

%% set directories and load data

basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
load(fullfile(basedir,['LSD_ls463_cat.mat'])) % make sure this file matches the split you want to run below
addpath(genpath('code'))


savedir = fullfile(basedir,'repkmeans');
mkdir(savedir); cd(savedir);


%% set inputs

concTS = TS_ls463_nomean(6294:end,:); %make sure this matches the split

split = 30; % i am using split to denote different processing applied to data 
% split 0 corresponds to gsr only (LSDgsr_cat.mat)
% split 1 corresponds to gsr + bp 0.008-0.09 (LSDgsr_bp_cat.mat)
% split 2 corresponds to gsr + bp 0.008-0.2 (LSDgsr_bp2_cat.mat) - 11/17
% accidentally overwrote k2,3 for split 2 with split 5 data!!
% split 3 corresponds to LSD ONLY gsr + bp 0.008-0.2  (LSDgsr_bp2_LSD.mat)
% split 4 corresponds to PL ONLY gsr + bp 0.008-0.2  (LSDgsr_bp2_PL.mat)
% split 5 corresponds to the data fully preprocessed by ICL, parcellated
% into sch200 +32Tian, and GSR applied after parcellation
% (LSD_sch232_cat.mat TS_gsr_sch232) 11/17 - now converted to double so
% kmeans should work better
% split 6 is same as 5 but only the LSD data (no PL)
% (TS_gsr_sch232(1:6293,:))
% split 7 is same as 5 but only the PL data (no LSD)
% (TS_gsr_sch232(6294:12586,:))
% split 8 is trying to run split 5 data while initializing the kmeans on
% centroids from split 7.
% split 9 is same as 5 but the first 90 regions of the aal116 atlas
% see ami_calc.m for full list of split descriptions

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
