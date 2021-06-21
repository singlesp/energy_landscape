%% script to get violin plot data - can alternatively use Cornblath functions
clear all; close all;clc
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
addpath(genpath('code'))
%% set inputs

split=22; % i am using split to denote different processing applied to data 
% see ami_calc.m for full list of split descriptions
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory

TR = 217;
nscans = 29; %only different from nsubjs if you have multiple scans per subj/conditions
LSD_stop = nscans;
PL_start = LSD_stop+1;
tot = nscans*2;
nsubjs=15;

scan_length = 7.33; %length of scan in minutes to get appearances per minute

subjInd=[repelem([1 3:nsubjs],1),repelem(1:nsubjs,1)]; % index data from each subject
% subjInd=[repelem(1:nsubjs,1)]; %for data with only 1 scan per subj/cond

%% big loop
for numClusters=[5]
    load(fullfile(savedir,['Partition_bp',num2str(split),'_k',num2str(numClusters),'.mat']))
    %% count number of each cluster per scan
    A=reshape(partition,TR,[]);
    count = zeros(numClusters,tot);
    for b=1:tot
        [count(:,b),~] = hist(A(:,b),1:numClusters);
    end
    
    %% Calculate Fractional Occupancy
    LSD=count(:,1:LSD_stop); 
    LSDfo1 = LSD/TR; %data for all scans
    
    LSDfo=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        LSDfo(:,i) = mean(LSDfo1(:,subjInd==i),2);
    end
    
    PL=count(:,PL_start:tot);
    PLfo1 = PL/TR;
    
    PLfo=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        PLfo(:,i) = mean(PLfo1(:,subjInd==i),2);
    end
    
    %% Calculate Dwell Time and Appearance Rate
    dwell = []; %zeros(numClusters,220,89);
    sumApp=0;
    for c=1:numClusters
        for b=1:tot
            appear=find(A(:,b)==c);
            [maxIndex,~]=size(appear);
            appear=vertcat(appear,zeros(1,1)); %so that while statement will allow us to count the last index
            s=1;
            i=1;
            a=1;
            sumApp=sumApp+maxIndex;
            while a<maxIndex+1
                if appear(a+1,1)==appear(a,1)+1 %if the next index value follows the current one in numerical order, then we are still dwelling in this state
                    s=s+1;
                    dwell(c,i,b)=s;
                    a=a+1;
                else
                    dwell(c,i,b)=s;
                    i=i+1;
                    a=a+1;
                    s=1;
                end
            end
            if sum(dwell(c,:,b)) ~= maxIndex
                disp(['Warning! Cluster ',num2str(c),' and Scan ',num2str(b),' sum does not match appearance count']);
            end
        end
    end
    
    
    LSD = dwell(:,:,1:LSD_stop); 
    LSD_count=zeros(numClusters,3,nscans); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    LSDdt1=[];
    LSDar1=[];
    [~,maxIndex,~]=size(dwell);
    for c=1:numClusters
        for b=1:nscans
            LSD_count(c,1,b)=LSD_count(c,1,b)+sum(LSD(c,:,b));
            a=1;
            while a<=maxIndex
                if LSD(c,a,b)~= 0
                    LSD_count(c,2,b)=LSD_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            LSDdt1(c,b)=LSD_count(c,1,b)/LSD_count(c,2,b)*2; %total time/#appear *2 to convert to seconds
            LSDar1(c,b)=LSD_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
    
    LSDdt=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        LSDdt(:,i) = mean(LSDdt1(:,subjInd==i),2);
    end
    
    LSDar=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        LSDar(:,i) = mean(LSDar1(:,subjInd==i),2);
    end
    
    
    
    PL=dwell(:,:,PL_start:tot);
    PL_count=zeros(numClusters,3,nscans); %rows represent each cluster,column 1=total time, 2=total number appearances,3=avg dwell time
    PLdt1=[];
    PLar1=[];
    [~,maxIndex,~]=size(dwell);
    for c=1:numClusters
        for b=1:nscans
            PL_count(c,1,b)=PL_count(c,1,b)+sum(PL(c,:,b));
            a=1;
            while a<=maxIndex
                if PL(c,a,b)~= 0
                    PL_count(c,2,b)=PL_count(c,2,b)+1;
                    a=a+1;
                else
                    a=a+1;
                end
            end
            PLdt1(c,b)=PL_count(c,1,b)/PL_count(c,2,b)*2; %total time/#appear *2 to convert to seconds
            PLar1(c,b)=PL_count(c,2,b)/scan_length; %appearance rate per minute = tot. appear / 7 min 20 s scan
        end
    end
    
    PLdt=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        PLdt(:,i) = mean(PLdt1(:,subjInd==i),2);
    end
    
    PLar=NaN(numClusters,nsubjs);
    
    %average scans by subject
    for i=1:nsubjs
        PLar(:,i) = mean(PLar1(:,subjInd==i),2);
    end
    
    %compute DT's and AR's averaged over clusters (does LSD overall have
    %lower DT's and higher AR's?)
    DT=NaN(1,nsubjs*2);
    DT(1:nsubjs)=mean(LSDdt,1);
    DT(nsubjs+1:nsubjs*2)=mean(PLdt,1);
    
%     [~,pDT,~,tDT] = ttest(DT(1:nsubjs),DT(nsubjs+1:nsubjs*2));
    
    AR=NaN(1,nsubjs*2);
    AR(1:nsubjs)=mean(LSDar,1);
    AR(nsubjs+1:nsubjs*2)=mean(PLar,1);
    
%     [~,pAR,~,tAR] = ttest(AR(1:nsubjs),AR(nsubjs+1:nsubjs*2));
    
    
    %% F-test for equality of variances
    
    %combine pre+post conditions into one group and during music into another
    
    % dvar = zeros(numClusters,3,2); %get a quick measure of %differences in variance for the two groups
    %
    % for i=1:numClusters
    %     dvar(i,1,1) = (var(LSDfo(i,[1:14 30:44]))-var(PLfo(i,[1:14 30:44])))/(var(LSDfo(i,[1:14 30:44])))*100;
    %     dvar(i,2,1) = (var(LSDdt(i,[1:14 30:44]))-var(PLdt(i,[1:14 30:44])))/(var(LSDdt(i,[1:14 30:44])))*100;
    %     dvar(i,3,1) = (var(LSDar(i,[1:14 30:44]))-var(PLar(i,[1:14 30:44])))/(var(LSDar(i,[1:14 30:44])))*100;
    %
    %     dvar(i,1,2) = (var(LSDfo(i,15:29))-var(PLfo(i,15:29)))/(var(LSDfo(i,15:29)))*100;
    %     dvar(i,2,2) = (var(LSDdt(i,15:29))-var(PLdt(i,15:29)))/(var(LSDdt(i,15:29)))*100;
    %     dvar(i,3,2) = (var(LSDar(i,15:29))-var(PLar(i,15:29)))/(var(LSDar(i,15:29)))*100;
    % end
    %
    % pvaluesPP = zeros(6,3);
    % pvaluesdur = zeros(6,3);
    
    %
    
    % for i=1:numClusters
    %     [~,pvaluesPP(i,1)] = vartest2(LSDfo(i,[1:14 30:44]),PLfo(i,[1:14 30:44]));
    %     [~,pvaluesPP(i,2)] = vartest2(LSDdt(i,[1:14 30:44]),PLdt(i,[1:14 30:44]));
    %     [~,pvaluesPP(i,3)] = vartest2(LSDar(i,[1:14 30:44]),PLar(i,[1:14 30:44]));
    %
    %     [~,pvaluesdur(i,1)] = vartest2(LSDfo(i,15:29),PLfo(i,15:29));
    %     [~,pvaluesdur(i,2)] = vartest2(LSDdt(i,15:29),PLdt(i,15:29));
    %     [~,pvaluesdur(i,3)] = vartest2(LSDar(i,15:29),PLar(i,15:29));
    % end
    %
    % pvaluesPP = reshape(mafdr(reshape(pvaluesPP,1,18),"BHFDR",1),6,3);
    % pvaluesdur = reshape(mafdr(reshape(pvaluesdur,1,18),"BHFDR",1),6,3);
    
    %% Detect outliers
    
    % outPP = NaN(6,29,6);
    % outdur = NaN(6,15,6);
    %
    % for i=1:numClusters
    %     outPP(i,:,1) = isoutlier(LSDfo(i,[1:14 30:44]));
    %     outPP(i,:,2) = isoutlier(LSDdt(i,[1:14 30:44]));
    %     outPP(i,:,3) = isoutlier(LSDar(i,[1:14 30:44]));
    %     outPP(i,:,4) = isoutlier(PLfo(i,[1:14 30:44]));
    %     outPP(i,:,5) = isoutlier(PLdt(i,[1:14 30:44]));
    %     outPP(i,:,6) = isoutlier(PLar(i,[1:14 30:44]));
    %
    %     outdur(i,:,1) = isoutlier(LSDfo(i,15:29));
    %     outdur(i,:,2) = isoutlier(LSDdt(i,15:29));
    %     outdur(i,:,3) = isoutlier(LSDar(i,15:29));
    %     outdur(i,:,4) = isoutlier(PLfo(i,15:29));
    %     outdur(i,:,5) = isoutlier(PLdt(i,15:29));
    %     outdur(i,:,6) = isoutlier(PLar(i,15:29));
    %
    % end
    %
    % sumPP(1,:) = sum(sum(outPP(:,:,1:3),3),1);
    % sumPP(2,:) = sum(sum(outPP(:,:,4:6),3),1);
    % sumdur(1,:) = sum(sum(outdur(:,:,1:3),3),1);
    % sumdur(2,:) = sum(sum(outdur(:,:,4:6),3),1);
    %
    
    %% Save
    
    clusters=char(clusterNames);
    save(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']), 'LSDfo', 'PLfo', 'LSDdt', 'PLdt', 'LSDar', 'PLar', 'clusters','DT','AR');%,'pDT','pAR','tDT','tAR')

end
