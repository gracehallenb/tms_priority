%% srt
par = parula
beta = 0.65
c1 = par(1,:,:)
c2 = par(2,:,:)
c3 = brighten(c1,beta)
c4 = brighten(c2,beta)
% c3 = par(7,:,:)  
% c4 = par(8,:,:)
c5 = par(17,:,:)
c6 = par(18,:,:)
c7 = brighten(c5, beta)
c8 = brighten(c6, beta)
% c7 = par(23,:,:)
% c8 = par(24,:,:)
c9 = par(34,:,:)
c10 = par(35,:,:)
c11 = brighten(c9,beta)
c12 = brighten(c10,beta)
% c11 = par(39,:,:)
% c12 = par(40,:,:)
c13 = par(55,:,:)
c14 = par(56,:,:)
c15 = brighten(c13,beta)
c16 = brighten(c14,beta)
% c15 = par(63,:,:)
% c16 = par(64,:,:)

parspec = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16]

%% pilot 
%ii_results.median_no_break_left_srt
%subj = {'subj01','subj02','subj03'};
%subj = {'subj01','subj02','subj03'};
subj = {'subj01'};
%cond  =  {'pilot','sham','l_spcs'}; %'l_ips2',
cond  =  {'pilot'};
num_cond = length(cond);
num_subj = length(subj);  
no_break_lo_left_srt_subj_pilot = [];
median_no_break_lo_left_srt_subj_pilot = [];
median_no_break_hi_left_srt_subj = [];
median_no_break_lo_right_srt_subj = [];
median_no_break_hi_right_srt_subj = [];

no_break_lo_right_srt_subj_pilot = [];
for ss = 1:length(subj);
    for cc= 1:length(cond);
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
        no_break_lo_left_srt_subj_pilot = [no_break_lo_left_srt_subj_pilot; resultsfile.ii_results_lo.no_break_left_srt]
        no_break_lo_right_srt_subj_pilot = [no_break_lo_right_srt_subj_pilot; resultsfile.ii_results_lo.no_break_right_srt];
         median_no_break_lo_left_srt_subj_pilot(ss,cc) =  resultsfile.ii_results_lo.median_no_break_left_srt;
         median_no_break_lo_right_srt_subj_pilot(ss,cc) =  resultsfile.ii_results_lo.median_no_break_right_srt;
    end
end

no_break_hi_left_srt_subj_pilot = []; 
no_break_hi_right_srt_subj_pilot = [];
for ss = 1:length(subj);
    for cc= 1:length(cond);
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
         no_break_hi_left_srt_subj_pilot = [no_break_hi_left_srt_subj_pilot; resultsfile.ii_results_hi.no_break_left_srt]
        no_break_hi_right_srt_subj_pilot = [no_break_hi_right_srt_subj_pilot; resultsfile.ii_results_hi.no_break_right_srt];
        
%         no_break_hi_left_srt_subj{ss,cc} = resultsfile.ii_results_hi.no_break_left_srt
%         no_break_hi_right_srt_subj{ss,cc} = resultsfile.ii_results_hi.no_break_right_srt
        median_no_break_hi_left_srt_subj_pilot(ss,cc) =  resultsfile.ii_results_hi.median_no_break_left_srt;
         median_no_break_hi_right_srt_subj_pilot(ss,cc) =  resultsfile.ii_results_hi.median_no_break_right_srt;
    end
end



for ii = 1:length(num_subj);
hileftsemsubj_pilot(:,ii) = std(no_break_hi_left_srt_subj_pilot(:,ii))/sqrt(length(no_break_hi_left_srt_subj_pilot(:,ii)))
hirightsemsubj_pilot(:,ii) = std(no_break_hi_right_srt_subj_pilot(:,ii))/sqrt(length(no_break_hi_right_srt_subj_pilot(:,ii)))
loleftsemsubj_pilot(:,ii) = std(no_break_lo_left_srt_subj_pilot(:,ii))/sqrt(length(no_break_lo_left_srt_subj_pilot(:,ii)))
lorightsemsubj_pilot(:,ii) = std(no_break_lo_right_srt_subj_pilot(:,ii))/sqrt(length(no_break_lo_right_srt_subj_pilot(:,ii)))
end 





for ii = 1:num_cond;
group_median_lo_right_pilot(ii) = median(median_no_break_lo_left_srt_subj_pilot(:,ii));
group_median_lo_left_pilot(ii) = median(median_no_break_lo_right_srt_subj_pilot(:,ii));
group_median_hi_left_pilot(ii) = median(median_no_break_hi_left_srt_subj_pilot(:,ii));
group_median_hi_right_pilot(ii) = median(median_no_break_hi_right_srt_subj_pilot(:,ii));
end

% for ii = 1:num_cond;
% hileftsem_pilot(:,ii) = std(median_no_break_hi_left_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_hi_left_srt_subj_pilot(:,ii)))
% hirightsem_pilot(:,ii) = std(median_no_break_hi_right_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_hi_right_srt_subj_pilot(:,ii)))
% loleftsem_pilot(:,ii) = std(median_no_break_lo_left_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_lo_left_srt_subj_pilot(:,ii)))
% lorightsem_pilot(:,ii) = std(median_no_break_lo_right_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_lo_right_srt_subj_pilot(:,ii)))
% end 

% %%
% p_change_lo_left_sham = abs(median_no_break_lo_left_srt_subj(:,1) - median_no_break_lo_left_srt_subj(:,2))/(median_no_break_lo_left_srt_subj(:,1)).*100
% p_change_lo_right_sham = abs(median_no_break_lo_right_srt_subj(:,1) - median_no_break_lo_right_srt_subj(:,2))/(median_no_break_lo_right_srt_subj(:,1)).*100
% p_change_hi_left_sham = abs(median_no_break_hi_left_srt_subj(:,1) - median_no_break_hi_left_srt_subj(:,2))/(median_no_break_hi_left_srt_subj(:,1)).*100
% p_change_hi_right_sham = abs(median_no_break_hi_right_srt_subj(:,1) - median_no_break_hi_right_srt_subj(:,2))/(median_no_break_hi_right_srt_subj(:,1)).*100
% 
% p_change_lo_left_lspcs = abs(median_no_break_lo_left_srt_subj(:,1) - median_no_break_lo_left_srt_subj(:,3))/(median_no_break_lo_left_srt_subj(:,1)).*100
% p_change_lo_right_lspcs = abs(median_no_break_lo_right_srt_subj(:,1) - median_no_break_lo_right_srt_subj(:,3))/(median_no_break_lo_right_srt_subj(:,1)).*100
% p_change_hi_left_lspcs = abs(median_no_break_hi_left_srt_subj(:,1) - median_no_break_hi_left_srt_subj(:,3))/(median_no_break_hi_left_srt_subj(:,1)).*100
% p_change_hi_right_lspcs = abs(median_no_break_hi_right_srt_subj(:,1) - median_no_break_hi_right_srt_subj(:,3))/(median_no_break_hi_right_srt_subj(:,1)).*100
% 
% % x = [median_no_break_hi_left_srt_subj(1,1) median_no_break_hi_left_srt_subj(1,1)]
% % y = [median_no_break_hi_left_srt_subj(1,2) median_no_break_hi_left_srt_subj(1,3)]
% % sz= 50;
% % figure;
% scatter(x,y, sz,'filled')

sem_vect_pilot = [hileftsemsubj_pilot hirightsemsubj_pilot loleftsemsubj_pilot lorightsemsubj_pilot]

pilot = [group_median_hi_left_pilot(1) group_median_hi_right_pilot(1) group_median_lo_left_pilot(1) group_median_lo_right_pilot(1)]; 
%sham = [group_median_hi_left(2) group_median_hi_right(2) group_median_lo_left(2) group_median_lo_right(2)];
%lips2 = [group_median_hi_left(3) group_median_hi_right(3) group_median_lo_left(3) group_median_lo_right(3)];
%lspcs = [group_median_hi_left(2) group_median_hi_right(2) group_median_lo_left(2) group_median_lo_right(2)];
all_srt_pilot = [pilot]; 
%%
subj = {'subj01'};
%cond  =  {'pilot','sham','l_spcs'}; %'l_ips2',
cond  =  {'l_spcs'};
num_cond = length(cond);
num_subj = length(subj);  
no_break_lo_left_srt_subj_spscs = [];
median_no_break_lo_left_srt_subj_spscs = [];
median_no_break_hi_left_srt_subj = [];
median_no_break_lo_right_srt_subj = [];
median_no_break_hi_right_srt_subj = [];

no_break_lo_right_srt_subj_spscs = [];
for ss = 1:length(subj);
    for cc= 1:length(cond);
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
        no_break_lo_left_srt_subj_spscs = [no_break_lo_left_srt_subj_spscs; resultsfile.ii_results_lo.no_break_left_srt]
        no_break_lo_right_srt_subj_spscs = [no_break_lo_right_srt_subj_spscs; resultsfile.ii_results_lo.no_break_right_srt];
         median_no_break_lo_left_srt_subj_spscs(ss,cc) =  resultsfile.ii_results_lo.median_no_break_left_srt;
         median_no_break_lo_right_srt_subj_spscs(ss,cc) =  resultsfile.ii_results_lo.median_no_break_right_srt;
    end
end

no_break_hi_left_srt_subj_spscs = []; 
no_break_hi_right_srt_subj_spscs = [];
for ss = 1:length(subj);
    for cc= 1:length(cond);
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
         no_break_hi_left_srt_subj_spscs = [no_break_hi_left_srt_subj_spscs; resultsfile.ii_results_hi.no_break_left_srt]
        no_break_hi_right_srt_subj_spscs = [no_break_hi_right_srt_subj_spscs; resultsfile.ii_results_hi.no_break_right_srt];
        
%         no_break_hi_left_srt_subj{ss,cc} = resultsfile.ii_results_hi.no_break_left_srt
%         no_break_hi_right_srt_subj{ss,cc} = resultsfile.ii_results_hi.no_break_right_srt
        median_no_break_hi_left_srt_subj_spscs(ss,cc) =  resultsfile.ii_results_hi.median_no_break_left_srt;
         median_no_break_hi_right_srt_subj_spscs(ss,cc) =  resultsfile.ii_results_hi.median_no_break_right_srt;
    end
end



for ii = 1:length(num_subj);
hileftsemsubj_spscs(:,ii) = std(no_break_hi_left_srt_subj_spscs(:,ii))/sqrt(length(no_break_hi_left_srt_subj_spscs(:,ii)))
hirightsemsubj_spscs(:,ii) = std(no_break_hi_right_srt_subj_spscs(:,ii))/sqrt(length(no_break_hi_right_srt_subj_spscs(:,ii)))
loleftsemsubj_spscs(:,ii) = std(no_break_lo_left_srt_subj_spscs(:,ii))/sqrt(length(no_break_lo_left_srt_subj_spscs(:,ii)))
lorightsemsubj_spscs(:,ii) = std(no_break_lo_right_srt_subj_spscs(:,ii))/sqrt(length(no_break_lo_right_srt_subj_spscs(:,ii)))
end 





for ii = 1:num_cond;
group_median_lo_right_spscs(ii) = median(median_no_break_lo_left_srt_subj_spscs(:,ii));
group_median_lo_left_spscs(ii) = median(median_no_break_lo_right_srt_subj_spscs(:,ii));
group_median_hi_left_spscs(ii) = median(median_no_break_hi_left_srt_subj_spscs(:,ii));
group_median_hi_right_spscs(ii) = median(median_no_break_hi_right_srt_subj_spscs(:,ii));
end

% for ii = 1:num_cond;
% hileftsem_spscs(:,ii) = std(median_no_break_hi_left_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_hi_left_srt_subj_pilot(:,ii)))
% hirightsem_pilot(:,ii) = std(median_no_break_hi_right_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_hi_right_srt_subj_pilot(:,ii)))
% loleftsem_pilot(:,ii) = std(median_no_break_lo_left_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_lo_left_srt_subj_pilot(:,ii)))
% lorightsem_pilot(:,ii) = std(median_no_break_lo_right_srt_subj_pilot(:,ii))/sqrt(length(median_no_break_lo_right_srt_subj_pilot(:,ii)))
% end 

% %%
% p_change_lo_left_sham = abs(median_no_break_lo_left_srt_subj(:,1) - median_no_break_lo_left_srt_subj(:,2))/(median_no_break_lo_left_srt_subj(:,1)).*100
% p_change_lo_right_sham = abs(median_no_break_lo_right_srt_subj(:,1) - median_no_break_lo_right_srt_subj(:,2))/(median_no_break_lo_right_srt_subj(:,1)).*100
% p_change_hi_left_sham = abs(median_no_break_hi_left_srt_subj(:,1) - median_no_break_hi_left_srt_subj(:,2))/(median_no_break_hi_left_srt_subj(:,1)).*100
% p_change_hi_right_sham = abs(median_no_break_hi_right_srt_subj(:,1) - median_no_break_hi_right_srt_subj(:,2))/(median_no_break_hi_right_srt_subj(:,1)).*100
% 
% p_change_lo_left_lspcs = abs(median_no_break_lo_left_srt_subj(:,1) - median_no_break_lo_left_srt_subj(:,3))/(median_no_break_lo_left_srt_subj(:,1)).*100
% p_change_lo_right_lspcs = abs(median_no_break_lo_right_srt_subj(:,1) - median_no_break_lo_right_srt_subj(:,3))/(median_no_break_lo_right_srt_subj(:,1)).*100
% p_change_hi_left_lspcs = abs(median_no_break_hi_left_srt_subj(:,1) - median_no_break_hi_left_srt_subj(:,3))/(median_no_break_hi_left_srt_subj(:,1)).*100
% p_change_hi_right_lspcs = abs(median_no_break_hi_right_srt_subj(:,1) - median_no_break_hi_right_srt_subj(:,3))/(median_no_break_hi_right_srt_subj(:,1)).*100
% 
% % x = [median_no_break_hi_left_srt_subj(1,1) median_no_break_hi_left_srt_subj(1,1)]
% % y = [median_no_break_hi_left_srt_subj(1,2) median_no_break_hi_left_srt_subj(1,3)]
% % sz= 50;
% % figure;
% scatter(x,y, sz,'filled')

sem_vect_spscs = [hileftsemsubj_spscs hirightsemsubj_spscs loleftsemsubj_spscs lorightsemsubj_spscs]

spscs = [group_median_hi_left_spscs(1) group_median_hi_right_spscs(1) group_median_lo_left_spscs(1) group_median_lo_right_spscs(1)]; 
%sham = [group_median_hi_left(2) group_median_hi_right(2) group_median_lo_left(2) group_median_lo_right(2)];
%lips2 = [group_median_hi_left(3) group_median_hi_right(3) group_median_lo_left(3) group_median_lo_right(3)];
%lspcs = [group_median_hi_left(2) group_median_hi_right(2) group_median_lo_left(2) group_median_lo_right(2)];
all_srt_spscs = [spscs]; 









%%
figure
subplot(1,2,1)
plot([1 2 3],median_no_break_hi_left_srt_subj(1,:),'o-','linewidth', 2) 
hold on; 
plot([1 2 3],median_no_break_hi_right_srt_subj(1,:),'o-','linewidth', 2) 
plot([1 2 3],median_no_break_lo_left_srt_subj(1,:),'o-','linewidth', 2) 
plot([1 2 3],median_no_break_lo_right_srt_subj(1,:),'o-','linewidth', 2) 
legend({'Hi L', 'Hi R', 'Lo L', 'Lo R'})
set(gca, 'Xtick', [1 2 3])
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 14)
ylim([100 300])
ylabel('Time (ms)')
title('Subject 1')
subplot(1,2,2)
plot([1 2 3],median_no_break_hi_left_srt_subj(2,:),'o-','linewidth', 2) 
hold on; 
plot([1 2 3],median_no_break_hi_right_srt_subj(2,:),'o-','linewidth', 2) 
plot([1 2 3],median_no_break_lo_left_srt_subj(2,:),'o-','linewidth', 2) 
plot([1 2 3],median_no_break_lo_right_srt_subj(2,:),'o-','linewidth', 2) 
legend({'Hi L', 'Hi R', 'Lo L', 'Lo R'})
set(gca, 'Xtick', [1 2 3])
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 14)
ylabel('Time (ms)')
title('Subject 2')
ylim([100 300])
% subplot(1,3,3)
% plot([1 2 3],median_no_break_hi_left_srt_subj(3,:),'o-','linewidth', 2) 
% hold on; 
% plot([1 2 3],median_no_break_hi_right_srt_subj(3,:),'o-','linewidth', 2) 
% plot([1 2 3],median_no_break_lo_left_srt_subj(3,:),'o-','linewidth', 2) 
% plot([1 2 3],median_no_break_lo_right_srt_subj(3,:),'o-','linewidth', 2) 
legend({'Hi L', 'Hi R', 'Lo L', 'Lo R'})
ylim([100 300])
xticklabel = {'Baseline', 'Sham', 'LsPCS'}
set(gca, 'Xtick', [1 2 3])
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 14)
ylabel('Time (ms)')
title('Subject 3')

%% redo plots by hi/lo on xaxis 



    figure(1);
    hold on;
    subplot(1,3,1)
    plot([1,2], [median_no_break_hi_left_srt_subj_pilot(ii,1)  median_no_break_lo_left_srt_subj_pilot(ii,1)],'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1) %ipsi
    hold on;
    plot([1,2], [median_no_break_hi_right_srt_subj_pilot(ii,1)  median_no_break_lo_right_srt_subj_pilot(ii,1)],'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1) %contra
%     plot([1,2], [median_no_break_hi_left_srt_subj(ii,2)  median_no_break_lo_left_srt_subj(ii,2)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13) %ipsi
%     hold on;
%     plot([1,2], [median_no_break_hi_right_srt_subj(ii,2)  median_no_break_lo_right_srt_subj(ii,2)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13)
    ylim([100 300])
    xlim([0 4])
    xti = {'High Priority', 'Low Priority'};
    set(gca, 'xtick', [1 2])
    set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
    legend('No TMS ipsi', 'No TMS contra', 'sPCS ipsi','sPCS contra')
    set(gca, 'fontsize', 14)
    errorbar(all_srt_pilot,sem_vect_pilot,'.', 'linewidth',1.5)
 figure(1);
    hold on;
    subplot(1,3,1)
    plot([1,2], [median_no_break_hi_left_srt_subj_spscs(ii,1)  median_no_break_lo_left_srt_subj_spscs(ii,1)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c1) %ipsi
    hold on;
    plot([1,2], [median_no_break_hi_right_srt_subj_spscs(ii,1)  median_no_break_lo_right_srt_subj_spscs(ii,1)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13) %contra
%     plot([1,2], [median_no_break_hi_left_srt_subj(ii,2)  median_no_break_lo_left_srt_subj(ii,2)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13) %ipsi
%     hold on;
%     plot([1,2], [median_no_break_hi_right_srt_subj(ii,2)  median_no_break_lo_right_srt_subj(ii,2)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13)
    ylim([100 300])
    xlim([0 4])
    xti = {'High Priority', 'Low Priority'};
    set(gca, 'xtick', [1 2])
    set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
    legend('No TMS ipsi', 'No TMS contra', 'sPCS ipsi','sPCS contra')
    set(gca, 'fontsize', 14)
    allsrt = [all_srt_pilot all_srt_spscs]
    allsem = [sem_vect_pilot sem_vect_spscs]
    errorbar(allsrt,allsem,'.', 'linewidth',1.5)
%%
subplot(1,3,2)
plot([1,2], [median_no_break_hi_left_srt_subj(2,1)  median_no_break_lo_left_srt_subj(2,1)],'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1) %ipsi
hold on;
plot([1,2], [median_no_break_hi_right_srt_subj(2,1)  median_no_break_lo_right_srt_subj(2,1)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13) %contra
ylim([100 300])
xlim([0 3])
xti = {'High Priority', 'Low Priority'}; 
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
legend('No TMS ipsi', 'No TMS contra', 'sPCS ipsi','sPCS contra')
set(gca, 'fontsize', 14)
hold on;
subplot(1,3,3)
plot([1,2], [median_no_break_hi_left_srt_subj(3,1)  median_no_break_lo_left_srt_subj(3,1)],'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1) %ipsi
hold on;
plot([1,2], [median_no_break_hi_right_srt_subj(3,1)  median_no_break_lo_right_srt_subj(3,1)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13) %contra
hold on;
ylim([100 300])
xlim([0 3])
xti = {'High Priority', 'Low Priority'}; 
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
legend('No TMS ipsi', 'No TMS contra', 'sPCS ipsi','sPCS contra')
set(gca, 'fontsize', 14)

% hold on;
% plot([1,2], [group_mean_hi_left(3)  group_mean_lo_left(3)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13) %ipsi
% hold on;

%%
figure
subplot(1,3,1)
plot([1 2],median_no_break_hi_left_srt_subj_pilot(1,:),'o-','linewidth', 2) 
hold on; 
plot([1 2],median_no_break_hi_right_srt_subj_pilot(1,:),'o-','linewidth', 2) 
plot([1 2],median_no_break_lo_left_srt_subj_pilot(1,:),'o-','linewidth', 2) 
plot([1 2],median_no_break_lo_right_srt_sub_pilotj(1,:),'o-','linewidth', 2) 
legend({'High ipsi', 'High contra', 'Low ipsi', 'Low contra'})
set(gca, 'Xtick', [1 2 3])
set(gca, 'Xticklabel',{'No TMS','LsPCS'}')
set (gca,'FontSize', 14)
ylim([100 300])
ylabel('Time (ms)')
title('Subject 1')
% subplot(1,3,2)
% plot([1 2],median_no_break_hi_left_srt_subj(2,:),'o-','linewidth', 2) 
% hold on; 
% plot([1 2],median_no_break_hi_right_srt_subj(2,:),'o-','linewidth', 2) 
% plot([1 2],median_no_break_lo_left_srt_subj(2,:),'o-','linewidth', 2) 
% plot([1 2],median_no_break_lo_right_srt_subj(2,:),'o-','linewidth', 2) 
% legend({'High ipsi', 'High contra', 'Low ipsi', 'Low contra'})
% set(gca, 'Xtick', [1 2 3])
% set(gca, 'Xticklabel',{'No TMS', 'LsPCS'}')
% set (gca,'FontSize', 14)
% ylabel('Time (ms)')
% title('Subject 2')
% ylim([100 300])
% subplot(1,3,3)
% plot([1 2],median_no_break_hi_left_srt_subj(3,:),'o-','linewidth', 2) 
% hold on; 
% plot([1 2],median_no_break_hi_right_srt_subj(3,:),'o-','linewidth', 2) 
% plot([1 2],median_no_break_lo_left_srt_subj(3,:),'o-','linewidth', 2) 
% plot([1 2],median_no_break_lo_right_srt_subj(3,:),'o-','linewidth', 2) 
% legend({'High ipsi', 'High contra', 'Low ipsi', 'Low contra'})
% ylim([100 300])
% xticklabel = {'Baseline','LsPCS'}
% set(gca, 'Xtick', [1 2 3])d
% set(gca, 'Xticklabel',{'No TMS', 'LsPCS'}')
% set (gca,'FontSize', 14)
% ylabel('Time (ms)')
% title('Subject 3')
% 
% 
% 
% %%
anova_group = [median_no_break_hi_left_srt_subj median_no_break_hi_right_srt_subj median_no_break_hi_left_srt_subj median_no_break_hi_right_srt_subj]
 anova_group = anova_group'
 anova_groupt = [anova_group(:,1); anova_group(:,2); anova_group(:,3)] 
% in this case, we flip the rows and columns to obtain each column as a
%subject
% stats setup
%subject1 = {'subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1'};
subject1 = {'subj1','subj1','subj1','subj1','subj1','subj1','subj1','subj1'};
subject2 = {'subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2'};
%subject2 = {'subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2','subj2'};
%subject3 = {'subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3'};
subject3 = {'subj3','subj3','subj3','subj3','subj3','subj3','subj3','subj3'};
subj = [subject1'; subject2'; subject3'];
%priority = {'high';'high';'low';'low';         'high';'high';'low';'low';       'high';'high';'low';'low'}
priority = {'high';'high';'low';'low'; 'high';'high';'low';'low'}
%hemi = {'left';'right';'left';'right';          'left';'right'; 'left'; 'right'; 'left'; 'right'; 'left'; 'right'}
hemi = {'left';'right';'left';'right';  'left'; 'right'; 'left'; 'right'}
% condition = {'pilot';'pilot';'pilot'; 'pilot'; 'sham' ;'sham'; 'sham'; 'sham'; 'lspcs';'lspcs';'lspcs';'lspcs'}
condition = {'pilot';'pilot';'pilot'; 'pilot'; 'lspcs';'lspcs';'lspcs';'lspcs'}
priority = [priority; priority; priority;]
hemi = [hemi;hemi;hemi];
condition = [condition; condition;condition];
[p tbl stats] = anovan(anova_groupt,{subj,priority,hemi,condition},'model','full','random',1,'varnames',{'subj','priority','hemi','condition'})

check = [subj priority hemi condition]






























figure(2)
bar(1,all_srt(1),'Facecolor',parspec(1,:,:),'EdgeColor',[1 1 1]);
    hold on;
bar(2,all_srt(2),'Facecolor',parspec(2,:,:),'EdgeColor',[1 1 1]);
bar(3,all_srt(3),'Facecolor',parspec(3,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(4,all_srt(4),'Facecolor',parspec(4,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(5,all_srt(5),'Facecolor',parspec(5,:,:),'EdgeColor',[1 1 1]);
bar(6,all_srt(6),'Facecolor',parspec(6,:,:),'EdgeColor',[1 1 1]);
bar(7,all_srt(7),'Facecolor',parspec(7,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(8,all_srt(8),'Facecolor',parspec(8,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(9,all_srt(9),'Facecolor',parspec(9,:,:),'EdgeColor',[1 1 1]);
bar(10,all_srt(10),'Facecolor',parspec(10,:,:),'EdgeColor',[1 1 1]);
bar(11,all_srt(11),'Facecolor',parspec(11,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(12,all_srt(12),'Facecolor',parspec(12,:,:),'EdgeColor',[1 1 1]); 
errorbar(all_srt,sem_vect,'.')

%% for 4 cond

figure
subplot(1,3,1)
plot([1 2 3 4],median_no_break_hi_left_srt_subj(1,:),'o-','linewidth', 2) 
hold on; 
plot([1 2 3 4],median_no_break_hi_right_srt_subj(1,:),'o-') 
plot([1 2 3 4],median_no_break_lo_left_srt_subj(1,:),'o-') 
plot([1 2 3 4],median_no_break_lo_right_srt_subj(1,:),'o-') 
legend({'hi left', 'hi right', 'lo left', 'lo right'})
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 12)
ylim([100 300])
subplot(1,3,2)
plot([1 2 3 4],median_no_break_hi_left_srt_subj(2,:),'o-') 
hold on; 
plot([1 2 3 4],median_no_break_hi_right_srt_subj(2,:),'o-') 
plot([1 2 3 4],median_no_break_lo_left_srt_subj(2,:),'o-') 
plot([1 2 3 4],median_no_break_lo_right_srt_subj(2,:),'o-') 
legend({'hi left', 'hi right', 'lo left', 'lo right'})
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 12)
ylim([100 300])
subplot(1,3,3)
plot([1 2 3 4],median_no_break_hi_left_srt_subj(3,:),'o-') 
hold on; 
plot([1 2 3 4],median_no_break_hi_right_srt_subj(3,:),'o-') 
plot([[1 2 3 4],median_no_break_lo_left_srt_subj(3,:),'o-') 
plot([1 2 3 4],median_no_break_lo_right_srt_subj(3,:),'o-') 
legend({'hi left', 'hi right', 'lo left', 'lo right'})
ylim([100 300])
all_srt = [pilot sham lspcs]; 
xticklabel = {'Baseline', 'Sham', 'LsPCS'}
set(gca, 'Xticklabel',{'No TMS', 'Sham', 'LsPCS'}')
set (gca,'FontSize', 12)

figure(2)
bar(1,all_srt(1),'Facecolor',parspec(1,:,:),'EdgeColor',[1 1 1]);
    hold on;
bar(2,all_srt(2),'Facecolor',parspec(2,:,:),'EdgeColor',[1 1 1]);
bar(3,all_srt(3),'Facecolor',parspec(3,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(4,all_srt(4),'Facecolor',parspec(4,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(5,all_srt(5),'Facecolor',parspec(5,:,:),'EdgeColor',[1 1 1]);
bar(6,all_srt(6),'Facecolor',parspec(6,:,:),'EdgeColor',[1 1 1]);
bar(7,all_srt(7),'Facecolor',parspec(7,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(8,all_srt(8),'Facecolor',parspec(8,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(9,all_srt(9),'Facecolor',parspec(9,:,:),'EdgeColor',[1 1 1]);
bar(10,all_srt(10),'Facecolor',parspec(10,:,:),'EdgeColor',[1 1 1]);
bar(11,all_srt(11),'Facecolor',parspec(11,:,:),'EdgeColor',[1 1 1]); %make lighter 
bar(12,all_srt(12),'Facecolor',parspec(12,:,:),'EdgeColor',[1 1 1]);