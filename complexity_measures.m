%% complexity measures (load a combined brain-state time-series)
clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs

split=22; % i am using split to denote different processing applied to data 
numClusters = 4;
savedir = fullfile(basedir,'results','example');mkdir(savedir);	

% load a combined brain-state time series (SOM+ and SOM- would be one
% state -- use comb_clusters.m to achieve this)
load(fullfile(savedir,['cPartition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'partition');
TR=217;
nscans=29;
tot = nscans*2;
nsubjs = 15;
subjInd=[repelem([1 3:nsubjs],1),repelem(1:nsubjs,1)]'; % index data from each subject

%% binarize partition

bin_partition = partition - 1; 

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
        LSDlz(:,i) = mean(LSDlz1(subjInd==i));
    end
    
    
    PLlz1 = LZ(nscans+1:tot);
    
    PLlz=NaN(1,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        PLlz(:,i) = mean(PLlz1(subjInd==i));
    end

[h,p,ci,t]=ttest(LSDlz,PLlz);

%% save

% save(fullfile(savedir,['LZcomplexity_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'PLlz','LSDlz');

%% shannon entropy of rss plots

% load sch454edgeFC.mat
% entropy = NaN(nscans,1);
% 
% for i=1:nscans
%     entropy(i) = wentropy(rss(:,i),'shannon');
% end
% 
% [h1,p1,ci1,t1]=ttest(entropy(1:14),entropy(30:43));
