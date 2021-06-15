%%%
% Make a yeo label list for the schaefer 200 + tian 32 atlas. 
% The Tian regions will be made an 8th sub-cortical network
%%%

%% 

clear all; close all;
parcDIR = '/Users/sps253/Documents/parcs/';
addpath(parcDIR);

%% create array of ROI string lables from text file

%load .txt file as table
atlaslabels = readtable([parcDIR,'labels/Schaefer2018_400Parcels_7Networks_order.txt']);

%extract the var2 which contains the label string
sch400 = atlaslabels(:,2);

%rename to var1 so can be combined with tian labels later
sch400.Properties.VariableNames={'Var1'};

%create blank table & fill with 'tian-subcortical' for the 32 subcort ROI's
tian54 = table('Size',[54 1],'VariableTypes',"string");
tian54(:,1)={'tian-subcortical'};

%creat blank table for the 232 rois and fill
ALarray = table('Size',[454 1],'VariableTypes',"string"); %'tian';
ALarray = [sch400;tian54];

%convert to array
ALarray = table2array(ALarray);

%% use regular expression matching to assign value 1-8 for the 7 yeo 
% networks plus an 8th subcortical for each of the 232 ROI's

nParc = 454;

network7labelsSch454=NaN(1,nParc);



for i = 1:nParc
    
    if regexpmatch(ALarray{i,1},['Vis']) %VIS
        network7labelsSch454(1,i)=1;
    elseif regexpmatch(ALarray{i,1},['SomMot']) %SOM
        network7labelsSch454(1,i)=2;
    elseif regexpmatch(ALarray{i,1},['DorsAttn']) %DAT
        network7labelsSch454(1,i)=3;
    elseif regexpmatch(ALarray{i,1},['SalVentAttn']) %VAT
        network7labelsSch454(1,i)=4;
    elseif regexpmatch(ALarray{i,1},['Limbic']) %LIM
        network7labelsSch454(1,i)=5;
    elseif regexpmatch(ALarray{i,1},['Cont']) %FPN
        network7labelsSch454(1,i)=6;
    elseif regexpmatch(ALarray{i,1},['Default']) %DMN
        network7labelsSch454(1,i)=7;
    elseif regexpmatch(ALarray{i,1},['tian']) %SUB
        network7labelsSch454(1,i)=8;
    end     
end


% csvwrite('/Users/sps253/Documents/ROI_maps/sch232_to_yeo.csv', network7labelsSch232);
csvwrite('/Users/sps253/Documents/ROI_maps/sch454_to_yeo.csv', network7labelsSch454);
