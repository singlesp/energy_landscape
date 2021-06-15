%%
% Start here with your analysis after you have chosen k: 
% load concatenated TS with [T, nparc] size.
% Specify numclusters
% 10 partitions are compared for mutual information
% the partition which scores the most mutual information with all 49 other
% cluseters is selected, plotted, and saved for further analysis

%requires AMI function added to path available here: https://www.mathworks.com/matlabcentral/fileexchange/33144-the-adjusted-mutual-information


%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000; % how many times to let kmeans try to converge before aborting rep
split = 22; % i am using split to denote different processing applied to data 


%% load BOLD data

% replace TS with your BOLD data formatted as a T-by-nparc matrix
% where T is the number of time points and nparc is the number of ROIs
load(fullfile(basedir,['LSD_ls463_cat.mat']),'TS_ls463_nomean') % make sure this file matches the split you want to run (see above)


TR=217;

concTS = TS_ls463_nomean; %make sure this matches the split

% concTS = zscore(concTS); %option to z-score

[T,nparc] = size(concTS);


%% generate N partitions that will be compared pair-wise for mutual information

for numClusters=[5 6] %if undecided, run loop. If confident about k, you can run over only 1 value for numClusters
    disp(['Starting clusters number: ',num2str(numClusters)]);
    N = 10; %number of partitions to create and compare
    parts = NaN(T,N); %partitions will be stored here
    D = NaN(N,T,numClusters); %distance matrices will be stored here
    
    %load centroids to initialize kmeans on (optional) 
    % MUST ADD <'Start', init> to the kmeans inputs to use this
%     load(['Partition_bp12_k',num2str(numClusters),'.mat'],'centroids'); %make sure file matches centroids you want to initialize from
%     init = NaN(numClusters,nparc,nreps);
%     for i=1:nreps
%         init(:,:,i) = centroids.';
%     end
    
    for i=1:N
        disp(['Clusters: ',num2str(numClusters),'. Starting kmeans number: ',num2str(i)]);
        [parts(:,i),~,~,D(i,:,:)] = kmeans(concTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxI);
    end
    
    %% calculate adjusted mutual information for every pair of partitions
    
    ami_results = NaN(N,N);
    
    for i=1:N
        for j=1:N
            ami_results(i,j) = ami(parts(:,i),parts(:,j));
        end
    end
    
    % assess
    [m,ind] = max(sum(ami_results,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
    partition = parts(:,ind); % take partition that has most agreement with all other for further analysis
    
    % plot
    f = figure;
    
    imagesc(ami_results); title(['Adjusted Mutal Information between Partitions k=',num2str(numClusters)]); colorbar;
    axis square; set(gca,'FontSize',8);
    f.PaperUnits = 'inches';
    f.PaperSize = [4 2];
    f.PaperPosition = [0 0 4 2];
    saveas(f,fullfile(savedir,['AMI_bp',num2str(split),'_k',num2str(numClusters),'.pdf']));
    
    
    %% compute centroids and plot
    
    centroids = GET_CENTROIDS(concTS,partition,numClusters);
    % name clusters based on alignment with Yeo resting state networks
%     centroids = zscore(centroids);
    clusterNames = NAME_CLUSTERS_ANGLE(centroids);  % need to add a prior Yeo partition labels for your parcellation
    [clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids);  % need to add a prior Yeo partition labels for your parcellation
    
    f = figure;
    subplot(1,2,1); imagesc(centroids); title('Centroids'); xticks(1:numClusters); xticklabels(clusterNames);
    colormap('plasma'); axis square; colorbar; set(gca,'FontSize',8); COLOR_TICK_LABELS(true,false,numClusters);
    subplot(1,2,2); imagesc(corr(centroids)); title('Centroid Similarity'); colorbar; caxis([-1 1]);
    colormap('plasma'); axis square; set(gca,'FontSize',8); xticks(1:numClusters); yticks(1:numClusters);
    xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
    COLOR_TICK_LABELS(true,true,numClusters);
    f.PaperUnits = 'inches';
    f.PaperSize = [4 2];
    f.PaperPosition = [0 0 4 2];
    saveas(f,fullfile(savedir,['Centroids_bp',num2str(split),'_k',num2str(numClusters),'.pdf']));
    
    %% save
    
    save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids','D');
end
