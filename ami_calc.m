%% step 1 after choosing k (or a range of k)

%%
clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);

split = 'main'

load(fullfile(['data/',split,'.mat']))

%% set inputs

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 50;	% how many times to repeat clustering. will choose lowest error solution
maxI = 1000; % how many times to let kmeans try to converge before aborting rep


[T,nparc] = size(concTS);


%% generate N partitions that will be compared pair-wise for mutual information

for numClusters=[4] %if undecided, run loop. If confident about k, you can run over only 1 value for numClusters
    disp(['Starting clusters number: ',num2str(numClusters)]);
    N = 10; %number of partitions to create and compare
    parts = NaN(T,N); %partitions will be stored here
    D = NaN(N,T,numClusters); %distance matrices will be stored here
    
    %load centroids from main analysis to initialize kmeans on (optional - useful for smaller data sets/replications such as psilo) 
    % MUST ADD <'Start', init> to the kmeans inputs to use this
%     load(['Partition_bpmain_k',num2str(numClusters),'.mat'],'centroids'); %make sure file matches centroids you want to initialize from
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
    
    % plot SI Figure 1
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
    
    %SI Figure 3
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
    
    save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids');
end
