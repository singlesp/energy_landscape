function [clusterNames,reorderClusters,clusterNamesSort,net7angle] = NAME_CLUSTERS_ANGLE(centroids)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

[nparc,numClusters] = size(centroids);

% if nparc > 400
%     load('data/yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
if nparc == 90 %|| nparc == 116 %will not work bc of sch116
    load('data/aal_to_yeo.csv'); 
    network7labels=aal_to_yeo;
    %network7labels(isnan(network7labels))=8; %replace NaN values with 8 so it is unassigned
    %network7labels=[ones(13,1); 2*ones(13,1); 3*ones(13,1); 4*ones(13,1); 5*ones(13,1); 6*ones(13,1); 7*ones(12,1)];
elseif nparc == 200 || nparc == 232
    load('data/sch232_to_yeo.csv'); 
    network7labels=sch232_to_yeo;
elseif nparc == 100 || nparc == 116
    load('data/sch116_to_yeo.csv'); 
    network7labels=sch116_to_yeo;
elseif nparc == 400 || nparc == 454
    load('data/sch454_to_yeo.csv');
    network7labels=sch454_to_yeo;
elseif nparc == 462
    load('data/Lausanne_463_subnetworks.mat');

    network7labels=subnetworks_reorder;
    network7labels(463)=[];
elseif nparc == 461 %ls463 with brainstem (last region) and region 14 removed (artefacts)
    load('data/Lausanne_463_subnetworks.mat');
    network7labels=subnetworks_reorder;
    network7labels([14 463])=[];
end

network7labels=network7labels(1:nparc);
network7labels=reshape(network7labels,nparc,1);

numNets = 7; %need to decide whether to include an 8th subcortical network
% make a matrix where each column corresponds to a labeled Yeo system in Lausanne parcellation
% the columns are binary vectors indicated whether a region belongs to corresponding Yeo system

binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);
% then duplicate this matrix, multiply by -1 and horizontally concatenate to
% provide separate names for when systems are low amplitude

binaryNetVectors = [binaryNetVectors, -1*binaryNetVectors];

YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

% calculate cosine of angle between binary state vector and centroids

net7angle = zeros(numClusters,numNets*2);

for K = 1:numClusters
    for B = 1:(numNets*2)
        net7angle(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    ind = find(net7angle(K,:) == max(net7angle(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesInit{K} = YeoNetNames{ind};
    if ind > numNets
        plusminus(K) = false;
    end
end

clusterNames = cellstr(clusterNamesInit);

%sort by name then plus-minus
%
[clusterNamesSort,I] = sort(clusterNamesInit);
[~,I2] = sort(plusminus(I));
clusterNamesSort = clusterNamesSort(I2);
reorderClusters = I(I2);
%}
