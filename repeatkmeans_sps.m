clear all; close all;clc

basedir = '/Users/sps253/Documents/energy_landscape'; %%<- your dir here!
cd(basedir);

%% load data

split = 'main' %%<- options: main, gsr, sch, music, psilo


load(fullfile(['data/',split,'.mat']))



savedir = fullfile(basedir,'repkmeans');
mkdir(savedir); cd(savedir);


%% set inputs

nreps = 50; % how many times to repeat clustering. will choose lowest error solution
distanceMethod = 'correlation';
maxI = 1000; % how many times you allow kmeans to try to converge

maxk = round(sqrt(TR)) -1; %k^2 must be less than TR to capture all transitions

for numClusters = 2:maxk %because we have 220 frames, k^2 must be less than 220 to capture all transitions
    
    disp(['K = ',num2str(numClusters)]);

    disp('start k-means');
    
    [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    save(fullfile(savedir,['kmeans',split,'k_',num2str(numClusters),'.mat']),'partition','sumd')
   
    clear partition
    clear sumd
    
    
end

disp('complete');
