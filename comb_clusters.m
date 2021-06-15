%% This script may be used to combine cluster identities in a partition, 
% also (more commonly) can be used to reorder and rename clusters in meta-state hierarchy
% there is a function in the Cornblath repo (reorderClusters.m) to do this automatically, however i did not investigate using it

clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs
numClusters = 5;
split=22; %must match file you want to load
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
TR = 217;
[X,~]=size(partition);

%% Reassign clusters in partition 

% example combination scenario
% {'SOM+';'DMN+';'SOM-';'DAT+';'DAT-';'DMN-'} <- important to look at what
% the cluster assignments are for a given partition before running this
% code. It will need to be adjusted based on what clusters you are
% combining and their partitcular assignemnts. 

% 1 and 3 both become 1 with new name SOM
% 2 and 6 become 2 with new name DMN
% 4 and 5 become 3 with new name DAT
% clusterNames = {'MS-1';'MS-2'}; % new combined cluster names

% Can also just reorder so that you have same cluster order across
% different analysis streams/splits of data

clusterNames = {'SOM-';'SOM+';'FPN-';'FPN+';}; %new cluster names in the
% new order

% numClusters = 4; %new number of clusters (combined)

c_partition = NaN(X,1); %initialize new partition (combined part.)

c_partition(partition == 1) = 1; %assign 1's to 1 
c_partition(partition == 3) = 2; %reassign 3's to 2 

c_partition(partition == 4) = 3; %assign 4's to 3 
c_partition(partition == 2) = 4; %reassign 2's to 4 

% if you are using higher k:
% c_partition(partition == 5) = 1; %assign 3's to 3 
% c_partition(partition == 3) = 6; %reassign 4's to 3 
% 
% c_partition(partition == 2) = 7; %assign 3's to 3 (LIM- to LIM)
% c_partition(partition == 4) = 8; %reassign 4's to 3 (LIM+ to LIM)

partition = c_partition; %rename to original naming for streamlined analysis.

%% reorder centroids (skip this section if combining instead of re-ordering)
[nparc,~]=size(centroids);
r_centroids = NaN(nparc,numClusters);

% Check to make sure the new order here matches the one above! (LH / RH
% sides of assignment equations should swap)
r_centroids(:,1) = centroids(:,1);
r_centroids(:,2) = centroids(:,3);
r_centroids(:,3) = centroids(:,4);
r_centroids(:,4) = centroids(:,2);

% % r_centroids(:,5) = centroids(:,2);
% % r_centroids(:,6) = centroids(:,3);
% % r_centroids(:,7) = centroids(:,2);
% % r_centroids(:,8) = centroids(:,4);

centroids = r_centroids;

%% save partition in separate file if combining (add prefix 'c' to filename) or overwrite old file if just re-ordering


save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids','D');
