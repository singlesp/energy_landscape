%% correlate energy calcs with other measures Figure 6

clear all; close all;
basedir = '/Users/sps253/Documents/energy_landscape';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory


%%
split='main' 
load(fullfile(['data/',split,'.mat']))
numClusters = 4;

T=0.001;


% load ratingsmeasures.mat FDdiff
% 
% if nsubjs==12
%     FDdiff([3 12 15]) = [];
% end

load(fullfile(savedir,['subjenergy_split',num2str(split),'_k',num2str(numClusters),'_T',num2str(T),'.mat']));

load(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']));


load(fullfile(savedir,['LZcomplexity_bp',num2str(split),'_k',num2str(numClusters),'.mat']));


%% Whole matrix mean diff and stdev diff

E_L = E_full(1:nsubjs,:);
E_P = E_full(nsubjs+1:nsubjs*2,:);


for i=1:nsubjs
    E_Li(i,:)=mean(reshape(E_L(i,:),[numClusters numClusters])');
    E_Pi(i,:)=mean(reshape(E_P(i,:),[numClusters numClusters])');
end

E_Di = (E_Li - E_Pi)./(E_Li+E_Pi);
E_Dirsp = reshape(E_Di,nsubjs*numClusters,1);

FO_diff = ((LSDfo - PLfo)./(LSDfo+PLfo))';
FO_diff_rsp = reshape(FO_diff,nsubjs*numClusters,1);

DT_diff = ((LSDdt - PLdt)./(LSDdt+PLdt))';
DT_diff_rsp = reshape(DT_diff,nsubjs*numClusters,1);

AR_diff = ((LSDar - PLar)./(LSDar+PLar))';
AR_diff_rsp = reshape(AR_diff,nsubjs*numClusters,1);

[r,p] = corr(E_Dirsp,FO_diff_rsp);
[r1,p1] = corr(E_Dirsp,DT_diff_rsp);
[r2,p2] = corr(E_Dirsp,AR_diff_rsp);

[r3,p3]=corr(E_Di,FO_diff);
[r4,p4]=corr(E_Di,DT_diff);
[r5,p5]=corr(E_Di,AR_diff);

%% partial corr with framewise displacement



E_L = E_full(1:nsubjs,:);
E_P = E_full(nsubjs+1:nsubjs*2,:);

E_L_avg = mean(E_L');
E_P_avg = mean(E_P');

E_L_std = std(E_L');
E_P_std = std(E_P');

diff_E_avg = (E_L_avg - E_P_avg) ./ (E_L_avg + E_P_avg);
diff_E_std = (E_L_std - E_P_std) ./ (E_L_std + E_L_std);

DT_L = mean(LSDdt);
DT_P = mean(PLdt);
AR_L = mean(LSDar);
AR_P = mean(PLar);
FO_L = mean(LSDfo);
FO_P = mean(PLfo);

diffDT = (DT_L - DT_P) ./ (DT_L+DT_P);
diffAR = (AR_L - AR_P) ./ (AR_L+AR_P);
diffFO = (FO_L - FO_P) ./ (FO_L+FO_P);

diffLZ = (LSDlz - PLlz) ./ (LSDlz+PLlz);

[r6,p6] = corr(diff_E_avg',diffFO');
[r7,p7] = corr(diff_E_avg',diffDT');
[r8,p8] = corr(diff_E_avg',diffAR');
[r12,p12] = corr(diff_E_avg',diffLZ');

[r9,p9] = corr(diff_E_std',diffFO');
[r10,p10] = corr(diff_E_std',diffDT');
[r11,p11] = corr(diff_E_std',diffAR');
[r13,p13] = corr(diff_E_std',diffLZ');

% [r6,p6] = partialcorr(diff_E_avg',diffFO',FDdiff);
% [r7,p7] = partialcorr(diff_E_avg',diffDT',FDdiff);
% [r8,p8] = partialcorr(diff_E_avg',diffAR',FDdiff);
% [r12,p12] = partialcorr(diff_E_avg',diffLZ',FDdiff);
% 
% [r9,p9] = partialcorr(diff_E_std',diffFO',FDdiff);
% [r10,p10] = partialcorr(diff_E_std',diffDT',FDdiff);
% [r11,p11] = partialcorr(diff_E_std',diffAR',FDdiff);
% [r13,p13] = partialcorr(diff_E_std',diffLZ',FDdiff);


%% plot Eavg and Estd vs DT, AR, and LZ

load colors

% myvar = repelem(1:4,15); % replace your color var for this
j=1;
f=figure;
for i=1:2
    if i==1
        my_y_vals = diff_E_avg;
        ylab='∆E_a_v_g';
    else
        my_y_vals = diff_E_std;
        ylab='∆E_s_t_d';
    end
    for k =1:3
        
        if k==1
            my_x_vals = diffDT;
            xlab='∆DT';
            if i==1
                rr=r7; pp=p7;
            else
                rr=r10; pp=p10;
            end
        elseif k==2
            my_x_vals = diffAR;
            xlab='∆AR';
            if i==1
                rr=r8; pp=p8;
            else
                rr=r11; pp=p11;
            end
        else 
            my_x_vals = diffLZ;
            xlab='∆LZ';
            if i==1
                rr=r12; pp=p12;
            else
                rr=r13; pp=p13;
            end
        end
        
        subplot(2,3,j)
        
        mksize = 50;
        scatter(my_x_vals, my_y_vals, mksize, 'filled')
        lsline
        ylabel(ylab); xlabel(xlab); title('LSD - PL');
        text(.02,-.06,['r = ',num2str(rr),' p = ' num2str(pp) ],'Color','b','Fontsize', 12);
        set(gca,'FontSize',18);
        set(gca,'TickLength',[0 0]);
        set(gca,'Fontname','arial');
        
        j=j+1;

    end
end



   
%% corr FD for each cond separately (SI Figures 20-22)

% %LSD
% [FOr,FOp] = corr(LSDfo',LSD_FD);
% [DTr,DTp] = corr(LSDdt',LSD_FD);
% [ARr,ARp]=corr(LSDar',LSD_FD);
% [Er,Ep]=corr(E_L_avg',LSD_FD);
% 
% figure;
% title('LSD');
% for i=1:numClusters
%    subplot(3,4,i)
%    scatter(LSDfo(i,:)',LSD_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['FO Cluster ',num2str(i)]);
%    text(0.25,0.1,{['R = ',num2str(round(FOr(i),4))], ['p = ',num2str(round(FOp(i),4))]});
%    
%    subplot(3,4,i+4)
%    scatter(LSDdt(i,:)',LSD_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['DT Cluster ',num2str(i)]);
%    text(7.5,0.1,{['R = ',num2str(round(DTr(i),4))], ['p = ',num2str(round(DTp(i),4))]});
%    
%    subplot(3,4,i+8)
%    scatter(LSDar(i,:)',LSD_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['AR Cluster ',num2str(i)]);
%    text(2.25,0.1,{['R = ',num2str(round(ARr(i),4))], ['p = ',num2str(round(ARp(i),4))]});
% end
% 
% figure;
% scatter(E_L_avg',LSD_FD,50,'filled'); lsline
% ylabel('FD'); xlabel(['E_a_v_g']);
%    text(1000000,0.1,{['R = ',num2str(round(Er,4))], ['p = ',num2str(round(Ep,4))]});
%    
% %%
% %PL
% 
% [FOr,FOp] = corr(PLfo',PL_FD);
% [DTr,DTp] = corr(PLdt',PL_FD);
% [ARr,ARp]=corr(PLar',PL_FD);
% [Er,Ep]=corr(E_P_avg',PL_FD);
% 
% figure;
% title('PL');
% for i=1:numClusters
%    subplot(3,4,i)
%    scatter(PLfo(i,:)',PL_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['FO Cluster ',num2str(i)]);
%    text(0.25,0.1,{['R = ',num2str(round(FOr(i),4))], ['p = ',num2str(round(FOp(i),4))]});
%    
%    subplot(3,4,i+4)
%    scatter(PLdt(i,:)',PL_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['DT Cluster ',num2str(i)]);
%    text(7.5,0.1,{['R = ',num2str(round(DTr(i),4))], ['p = ',num2str(round(DTp(i),4))]});
%    
%    subplot(3,4,i+8)
%    scatter(PLar(i,:)',PL_FD,50,'filled'); lsline
%    ylabel('FD'); xlabel(['AR Cluster ',num2str(i)]);
%    text(2.25,0.1,{['R = ',num2str(round(ARr(i),4))], ['p = ',num2str(round(ARp(i),4))]});
% end
% 
% figure;
% scatter(E_P_avg',PL_FD,50,'filled'); lsline
% ylabel('FD'); xlabel(['E_a_v_g']);
%    text(1000000,0.1,{['R = ',num2str(round(Er,4))], ['p = ',num2str(round(Ep,4))]});

