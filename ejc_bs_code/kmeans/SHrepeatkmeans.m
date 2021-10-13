clear all;close all;
split=12;
rng(split);
nreps = 500;
addpaths;
% load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
% load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
basedir = '/Users/sps253/Documents/brain_states-master';

load(fullfile(basedir,['LSD_sch454_cat.mat']),'TS_sch454_nomean') % make sure this file matches the split you want to run below
addpath(genpath('code'))


savedir = fullfile(basedir,'repkmeans','splithalves');
mkdir(savedir); cd(savedir);

%csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
concTS=TS_sch454_nomean;
% if zdim > 0
%     concTS = zscore(concTS,[],zdim);
% end
[nobs,nparc]=size(concTS);
disp('data loaded');
numClusters=4;
nsubjs=15;
TR=217;
subjInd=[repelem([1 3:nsubjs],TR),repelem(1:nsubjs,TR),repelem([1 3:nsubjs],TR),repelem(1:nsubjs,TR)]'; % index data from each subject

disp(['K = ',num2str(numClusters),'Split = ',num2str(split)]);
distanceMethod='correlation';
% savedir = ['SHkmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
% mkdir(savedir);
%%
cd(savedir); 

disp('start k-means');
for R = 1:nreps
	obs1 = randperm(nobs,floor(0.5*nobs));
    [partition1,~,sumd1] = kmeans(concTS(obs1,:),numClusters,'Distance',distanceMethod);
    partition1 = int8(partition1);

    obs2 = find(~ismember(1:nobs,obs1));	% cluster the other half
    [partition2,~,sumd2] = kmeans(concTS(obs2,:),numClusters,'Distance',distanceMethod);
    partition2 = int8(partition2); 
    save(['SHkmeans',num2str(split),'k_',num2str(numClusters),'rep',num2str(R),'.mat'],'partition1','sumd1','obs1','partition2','sumd2','obs2');
    clear partition1; clear partition2;
    clear sumd1; clear sumd2;
    disp(['K-means ',num2str(R)]);
end

disp('complete');
