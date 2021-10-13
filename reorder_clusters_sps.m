%% This script may be used to combine cluster identities in a partition, 
% also (more commonly) can be used to reorder and rename clusters in meta-state hierarchy
% there is a function in the Cornblath repo (reorderClusters.m) to do this automatically, however i did not investigate using it

% use max correlation to a set of centroids or visualize with
% systems_plot.sps first to determin new order and manually input order
% below

clear all; close all;clc
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);

%% set inputs
numClusters = 4;
split='main' %must match file you want to load

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))

[X,~]=size(partition);

%% Reassign clusters in partition 


clusterNames = {'MS-1a';'MS-1b';'MS-2a';'MS-2b';}; %new cluster names in the
% new order


c_partition = NaN(X,1); %initialize new partition 

c_partition(partition == 3) = 1; %assign 1's to 1 
c_partition(partition == 2) = 2; %reassign 3's to 2 

c_partition(partition == 1) = 3; %assign 4's to 3 
c_partition(partition == 4) = 4; %reassign 2's to 4 

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
r_centroids(:,1) = centroids(:,3);
r_centroids(:,2) = centroids(:,2);
r_centroids(:,3) = centroids(:,1);
r_centroids(:,4) = centroids(:,4);

% % r_centroids(:,5) = centroids(:,2);
% % r_centroids(:,6) = centroids(:,3);
% % r_centroids(:,7) = centroids(:,2);
% % r_centroids(:,8) = centroids(:,4);

centroids = r_centroids;

%% save partition in separate file if combining (add prefix 'c' to filename) or overwrite old file if just re-ordering


save(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']),'parts', 'ami_results', 'ind', 'partition', 'clusterNames', 'centroids');
