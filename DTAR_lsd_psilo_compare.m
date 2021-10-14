%% SI Figure 5
clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory


%% load LSD data
split ='main'
numClusters = 4;

load(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']));

%%

meanLSDdt = mean(LSDdt);
meanPLdt = mean(PLdt);
[h,p,~,t] = ttest(meanLSDdt,meanPLdt);

meanLSDar = mean(LSDar)
meanPLar = mean(PLar)
[h1,p1,~,t1] = ttest(meanLSDar,meanPLar);

figure; violin([meanLSDdt',meanPLdt']); ylim([4 13])

figure; violin([meanLSDar',meanPLar']); ylim([1 3])

%% load PSILO data
split = 'psilo'
numClusters = 4;


load(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']));

%%

meanLSDdt = mean(LSDdt);
meanPLdt = mean(PLdt);
[h,p,~,t] = ttest(meanLSDdt,meanPLdt);

meanLSDar = mean(LSDar)
meanPLar = mean(PLar)
[h1,p1,~,t1] = ttest(meanLSDar,meanPLar);

figure; violin([meanLSDdt',meanPLdt']); ylim([4 13])

figure; violin([meanLSDar',meanPLar']); ylim([1 3])