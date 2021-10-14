%% complexity measures (load a combined brain-state time-series) must have partition ordered in meta-state
% order i.e. (1 and 2 belong to MS-1; 3 and 4 belong to MS-2)
% the data generated here are used in E_corrs.m

clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);	
%% set inputs

split='psilo' 
load(fullfile(['data/',split,'.mat']))
numClusters = 4;


% load a combined brain-state time series (SOM+ and SOM- would be one
% state -- use comb_clusters.m to achieve this)
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'partition');

tot = nscans*2;


%% binarize partition

bin_partition = NaN(size(partition));
bin_partition(partition<3)=0;
bin_partition(partition>2)=1;


%% calc LZ complexity with LZ76 (function calc_lz_complexity.m from mathworks Copyright (c) 2012, Quang Thai)

start=1;
stop=TR;

LZ=NaN(tot,1);

for i=1:tot
    LZ(i) = calc_lz_complexity(bin_partition(start:stop),'exhaustive',1);
    start=start+TR;
    stop=stop+TR;
end

LSDlz1 = LZ(1:nscans);
LSDlz=NaN(1,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        LSDlz(:,i) = mean(LSDlz1(subj_scanInd==i));
    end
    
    
    PLlz1 = LZ(nscans+1:tot);
    
    PLlz=NaN(1,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        PLlz(:,i) = mean(PLlz1(subj_scanInd==i));
    end

[h,p,ci,t]=ttest(LSDlz,PLlz);

%% save

save(fullfile(savedir,['LZcomplexity_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'PLlz','LSDlz');


