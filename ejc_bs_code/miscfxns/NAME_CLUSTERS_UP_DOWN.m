function [clusterNamesUp,clusterNamesDown,net7angle_Up,net7angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

[nparc,numClusters] = size(centroids);

% if nparc > 400
%     load('data/yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
if nparc == 90 %|| nparc == 116 %116 will not work now with sch116 OK because don't use aal116  anyways
   load('/Users/sps253/Documents/ROI_maps/aal_to_yeo.csv'); 
   network7labels=aal_to_yeo;
   %network7labels(isnan(network7labels))=8; %replace NaN values with 8 so it is unassigned
   % network7labels=[ones(13,1); 2*ones(13,1); 3*ones(13,1); 4*ones(13,1); 5*ones(13,1); 6*ones(13,1); 7*ones(12,1)];
elseif nparc ==200 || nparc ==232
   load('/Users/sps253/Documents/ROI_maps/sch232_to_yeo.csv'); 
   network7labels=sch232_to_yeo;
elseif nparc == 100 || nparc == 116
    load('/Users/sps253/Documents/ROI_maps/sch116_to_yeo.csv'); 
    network7labels=sch116_to_yeo;
elseif nparc == 400 || nparc == 454
    load('/Users/sps253/Documents/ROI_maps/sch454_to_yeo.csv');
    network7labels=sch454_to_yeo;
elseif nparc == 462
    load('/Users/sps253/Documents/ROI_maps/Lausanne_463_subnetworks.mat');
    network7labels=subnetworks;
    network7labels=subnetworks_reorder;
    network7labels(463)=[];
elseif nparc == 461 %ls463 with brainstem (last region) and region 14 removed (artefacts)
    load('/Users/sps253/Documents/ROI_maps/Lausanne_463_subnetworks.mat');
    network7labels=subnetworks_reorder;
    network7labels([14 463])=[];
end

network7labels=network7labels(1:nparc);
network7labels=reshape(network7labels,nparc,1);

numNets = 7; %need to decide whether to include 8th network
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'}; %SUB = subcortical regions

% calculate cosine of angle between binary state vector and centroids

centroids_up = centroids .* (centroids > 0);
centroids_down = -1 * centroids .* (centroids < 0);     % make negative activity positive and get rid of positive activity

net7angle_Up = zeros(numClusters,numNets);
net7angle_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        net7angle_Up(K,B) = dot(centroids_up(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
        net7angle_Down(K,B) = dot(centroids_down(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesUp = cell(numClusters,1);
clusterNamesDown = cell(numClusters,1);
for K = 1:numClusters       % for up and down separately, calculate closest network
    Up_ind = find(net7angle_Up(K,:) == max(net7angle_Up(K,:)));
    Down_ind = find(net7angle_Down(K,:) == max(net7angle_Down(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesUp{K} = [YeoNetNames{Up_ind},'+'];
    clusterNamesDown{K} = [YeoNetNames{Down_ind},'-'];
end