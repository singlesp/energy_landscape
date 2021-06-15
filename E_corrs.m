%% correlate energy calcs with other measures

clear all; close all;
basedir = '/Users/sps253/Documents/brain_states-master';
cd(basedir);
savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory


%%
split = 22;
numClusters = 4;

load(fullfile(savedir,['subjcentroids_split',num2str(split),'_k',num2str(numClusters),'.mat']));

load(fullfile(savedir,['ViolinData_bp',num2str(split),'_k',num2str(numClusters),'.mat']));


load(fullfile(savedir,['LZcomplexity_bp',num2str(split),'_k',num2str(numClusters),'.mat']));


%% Whole matrix mean diff and stdev diff

E_L = E_full(1:15,:);
E_P = E_full(16:30,:);


for i=1:15
    E_Li(i,:)=mean(reshape(E_L(i,:),[numClusters numClusters])');
    E_Pi(i,:)=mean(reshape(E_P(i,:),[numClusters numClusters])');
end

E_Di = (E_Li - E_Pi)./(E_Li+E_Pi);
E_Dirsp = reshape(E_Di,60,1);

FO_diff = ((LSDfo - PLfo)./(LSDfo+PLfo))';
FO_diff_rsp = reshape(FO_diff,60,1);

DT_diff = ((LSDdt - PLdt)./(LSDdt+PLdt))';
DT_diff_rsp = reshape(DT_diff,60,1);

AR_diff = ((LSDar - PLar)./(LSDar+PLar))';
AR_diff_rsp = reshape(AR_diff,60,1);

[r,p] = corr(E_Dirsp,FO_diff_rsp);
[r1,p1] = corr(E_Dirsp,DT_diff_rsp);
[r2,p2] = corr(E_Dirsp,AR_diff_rsp);

[r3,p3]=corr(E_Di,FO_diff);
[r4,p4]=corr(E_Di,DT_diff);
[r5,p5]=corr(E_Di,AR_diff);

%% partial corr with framewise displacement

load ratingsmeasures.mat FDdiff

E_L = E_full(1:15,:);
E_P = E_full(16:30,:);

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

% [r6,p6] = corr(diff_E_avg',diffFO');
% [r7,p7] = corr(diff_E_avg',diffDT');
% [r8,p8] = corr(diff_E_avg',diffAR');
% 
% [r9,p9] = corr(diff_E_std',diffFO');
% [r10,p10] = corr(diff_E_std',diffDT');
% [r11,p11] = corr(diff_E_std',diffAR');

[r6,p6] = partialcorr(diff_E_avg',diffFO',FDdiff);
[r7,p7] = partialcorr(diff_E_avg',diffDT',FDdiff);
[r8,p8] = partialcorr(diff_E_avg',diffAR',FDdiff);
[r12,p12] = partialcorr(diff_E_avg',diffLZ',FDdiff);

[r9,p9] = partialcorr(diff_E_std',diffFO',FDdiff);
[r10,p10] = partialcorr(diff_E_std',diffDT',FDdiff);
[r11,p11] = partialcorr(diff_E_std',diffAR',FDdiff);
[r13,p13] = partialcorr(diff_E_std',diffLZ',FDdiff);


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

%% radial values of LIMBIC is MS-1a (high-amplitude)

LSD_lim_UP = net7angle_Up(1:15,1,5);
PL_lim_UP = net7angle_Up(16:30,1,5);

diffLIM_UP = (LSD_lim_UP - PL_lim_UP)./(LSD_lim_UP + PL_lim_UP);



%% Subjective measures

load ratingsmeasures.mat ratings RatingsTable DASCTABLE asc


[rr1,pp1] = partialcorr(diff_E_avg',ratings,FDdiff);
[rr2,pp2] = partialcorr(diff_E_std',ratings,FDdiff);
[rr3,pp3] = partialcorr(E_L_avg',ratings,FDdiff);
[rr4,pp4] = partialcorr(E_L_std',ratings,FDdiff);
[rr5,pp5] = partialcorr(E_P_avg',ratings,FDdiff);
[rr6,pp6] = partialcorr(E_P_std',ratings,FDdiff);

[rr7,pp7] = partialcorr(diffLZ',ratings,FDdiff);

[rr8,pp8] = partialcorr(diffLIM_UP,ratings,FDdiff);


[rrr1,ppp1] = partialcorr(diff_E_avg',asc,FDdiff);
[rrr2,ppp2] = partialcorr(diff_E_std',asc,FDdiff);
[rrr3,ppp3] = partialcorr(E_L_avg',asc,FDdiff);
[rrr4,ppp4] = partialcorr(E_L_std',asc,FDdiff);
[rrr5,ppp5] = partialcorr(E_P_avg',asc,FDdiff);
[rrr6,ppp6] = partialcorr(E_P_std',asc,FDdiff);

[rrr7,ppp7] = partialcorr(diffLZ',asc,FDdiff);

[rrr8,ppp8] = partialcorr(diffLIM_UP,asc,FDdiff);

% [rr2,pp2] = partialcorr(diff_E_std',ratings,FDdiff);
% [rr3,pp3] = partialcorr(E_L_avg',ratings,FDdiff);
% [rr4,pp4] = partialcorr(E_L_std',ratings,FDdiff);
% [rr5,pp5] = partialcorr(E_P_avg',ratings,FDdiff);
% [rr6,pp6] = partialcorr(E_P_std',ratings,FDdiff);

% [rr1,pp1] = corr(diff_E_avg',ratings);
% [rr2,pp2] = corr(diff_E_std',ratings);
% [rr3,pp3] = corr(E_L_avg',ratings);
% [rr4,pp4] = corr(E_L_std',ratings);

titles = {'Intensity', 'Complex Imagery','Simple Halls','Emotional Arousal','Positive Mood','Ego Dissolution'};

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(diff_E_avg',ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('∆E_a_v_g');
end

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(diffLZ',ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('∆LZ');
end

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(E_L_avg',ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('E\_L_a_v_g');
end

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(E_P_avg',ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('E\_P_a_v_g');
end

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(E_L_std',ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('E\_L_s_t_d');
end

figure;
for i=1:6
    subplot(3,2,i)
    
    scatter(diffLIM_UP,ratings(:,i),50,'filled');lsline
    title(titles{i});
    ylabel('Rating'); xlabel('∆LIM');
end

figure;
for i=1:11
    subplot(3,4,i)
    
    scatter(E_P_avg,asc(:,i),50,'filled');lsline
%     title(titles{i});
    ylabel('Rating'); xlabel('E_P_L');
end



%% trying to see if post-acute E reductions are present with un-paired comparisons of who got LSD on V1 vs V2

orderIND = [1 0 1 0 0 1 0 1 0 1 1 1 0 1 0]; %1 recieve LSD on V1, 0 otherwise

diffELSD = diff_E_avg(orderIND==1);
diffEPL=diff_E_avg(orderIND==0);
diffEPL(1,8)=NaN;


PLELSD = E_L_avg(orderIND==1);
PLEPL = E_P_avg(orderIND==0);
PLELSD = E_P_avg(orderIND==1);
PLEPL(1,8)=NaN;

figure; violin([diffELSD',diffEPL']); xticks([1 2]); xticklabels({'LSD v1','PL v1'}); title('Energy reduction');
figure; violin([PLELSD',PLEPL']); xticks([1 2]); xticklabels({'LSD v1','PL v1'}); title('PL Energy');
