%% SUBJ across COND --   PRIMARY EYE POSITION
%need to be able to access means of each subj close all
%% colors
par = parula;
beta = 0.65;
c1 = par(1,:,:);
c2 = par(2,:,:);
c3 = brighten(c1,beta);
c4 = brighten(c2,beta);

% c3 = par(7,:,:)  
% c4 = par(8,:,:)

c5 = par(17,:,:);
c6 = par(18,:,:);
c7 = brighten(c5, beta);
c8 = brighten(c6, beta);

% c7 = par(23,:,:)
% c8 = par(24,:,:)

c9 = par(34,:,:);
c10 = par(35,:,:);
c11 = brighten(c9,beta);
c12 = brighten(c10,beta);

% c11 = par(39,:,:);
% c12 = par(40,:,:);

c13 = par(55,:,:);
c14 = par(56,:,:);
c15 = brighten(c13,beta);
c16 = brighten(c14,beta);



% c15 = par(63,:,:)
% c16 = par(64,:,:)

parspec = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16];

%%
subj = {'subj01','subj02','subj03','subj04'}; %,'subj04'
cond  =  {'pilot','l_spcs','l_ips2'}; %'l_ips2',
num_cond = length(cond);
num_subj = length(subj);  

for ss = 1:num_subj;
    for cc = 1:num_cond;
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
        resultsfile = load(filename);
        median_primary_err_lo_left_cond(ss,cc) = resultsfile.ii_results_lo.median_no_break_left_primary_err_z_new;
        median_primary_err_lo_right_cond(ss,cc) = resultsfile.ii_results_lo.median_no_break_right_primary_err_z_new;
    end 
    for cc = 1:num_cond;
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
        median_primary_err_hi_left_cond(ss,cc) = resultsfile.ii_results_hi.median_no_break_left_primary_err_z_new;
        median_primary_err_hi_right_cond(ss,cc) = resultsfile.ii_results_hi.median_no_break_right_primary_err_z_new;
    end
end 

for ii = 1:num_cond;
group_mean_lo_right(ii) = mean(median_primary_err_lo_right_cond(:,ii)); %takes mean down columns (aka, each row is a subj, so mean of all subj in that cond)
group_mean_lo_left(ii) = mean(median_primary_err_lo_left_cond(:,ii));
group_mean_hi_left(ii) = mean(median_primary_err_hi_left_cond(:,ii));
group_mean_hi_right(ii) = mean(median_primary_err_hi_right_cond(:,ii));
end 

for ii = 1:num_cond;
all_cond = [group_mean_hi_left group_mean_hi_right group_mean_lo_left group_mean_lo_right];
end

for jj = 1:num_subj;
    subj_mean_pri_hi_left = mean(median_primary_err_hi_left_cond,2)
    subj_mean_pri_hi_right = mean(median_primary_err_hi_right_cond,2)
    subj_mean_pri_lo_left = mean(median_primary_err_lo_left_cond,2)
    subj_mean_pri_lo_right = mean(median_primary_err_lo_right_cond,2)
end


%% subject wise plotting 11/15 

subj01_demeand = mean([subj_mean_pri_hi_left(1) subj_mean_pri_hi_right(1) subj_mean_pri_lo_left(1) subj_mean_pri_lo_right(1)]);
subj02_demeand = mean([subj_mean_pri_hi_left(2) subj_mean_pri_hi_right(2) subj_mean_pri_lo_left(2) subj_mean_pri_lo_right(2)]); 
subj03_demeand = mean([subj_mean_pri_hi_left(3) subj_mean_pri_hi_right(3) subj_mean_pri_lo_left(3) subj_mean_pri_lo_right(3)]);
subj04_demeand = mean([subj_mean_pri_hi_left(4) subj_mean_pri_hi_right(4) subj_mean_pri_lo_left(4) subj_mean_pri_lo_right(4)]);

% cond = {'hi', 'lo'}
% hemi = {'right', 'left'}; 
% subj = {'subj01','subj02','subj03','subj04'};
% 
% for ii=1:length(subjects);
%     for jj = 1:length(cond);
%     subj_demeaned{ii,jj} = median_primary_err_hi_left_cond(ii,jj)- subj01_demeand
%     
%     
%     
% end

 subj01_phl_one_demeaned = median_primary_err_hi_left_cond(1,1)- subj01_demeand;
 subj01_phl_two_demeaned = median_primary_err_hi_left_cond(1,2)- subj01_demeand;
 subj01_phl_three_demeaned = median_primary_err_hi_left_cond(1,3)- subj01_demeand;
 subj01_phr_one_demeaned = median_primary_err_hi_right_cond(1,1)- subj01_demeand;
 subj01_phr_two_demeaned = median_primary_err_hi_right_cond(1,2)- subj01_demeand;
 subj01_phr_three_demeaned = median_primary_err_hi_right_cond(1,3)- subj01_demeand;

 subj01_pll_one_demeaned = median_primary_err_lo_left_cond(1,1)- subj01_demeand;
 subj01_pll_two_demeaned = median_primary_err_lo_left_cond(1,2)- subj01_demeand;
 subj01_pll_three_demeaned = median_primary_err_lo_left_cond(1,3)- subj01_demeand;
 subj01_plr_one_demeaned = median_primary_err_lo_right_cond(1,1)- subj01_demeand;
 subj01_plr_two_demeaned = median_primary_err_lo_right_cond(1,2)- subj01_demeand;
 subj01_plr_three_demeaned = median_primary_err_lo_right_cond(1,3)- subj01_demeand;
 
 
 subj02_phl_one_demeaned = median_primary_err_hi_left_cond(2,1)- subj02_demeand;
 subj02_phl_two_demeaned = median_primary_err_hi_left_cond(2,2)- subj02_demeand;
 subj02_phl_three_demeaned = median_primary_err_hi_left_cond(2,3)- subj02_demeand;
 subj02_phr_one_demeaned = median_primary_err_hi_right_cond(2,1)- subj02_demeand;
 subj02_phr_two_demeaned = median_primary_err_hi_right_cond(2,2)- subj02_demeand;
 subj02_phr_three_demeaned = median_primary_err_hi_right_cond(2,3)- subj02_demeand;

 subj02_pll_one_demeaned = median_primary_err_lo_left_cond(2,1)- subj02_demeand;
 subj02_pll_two_demeaned = median_primary_err_lo_left_cond(2,2)- subj02_demeand;
 subj02_pll_three_demeaned = median_primary_err_lo_left_cond(2,3)- subj02_demeand;
 subj02_plr_one_demeaned = median_primary_err_lo_right_cond(2,1)- subj02_demeand;
 subj02_plr_two_demeaned = median_primary_err_lo_right_cond(2,2)- subj02_demeand;
 subj02_plr_three_demeaned = median_primary_err_lo_right_cond(2,3)- subj02_demeand;
 
 subj03_phl_one_demeaned = median_primary_err_hi_left_cond(3,1)- subj03_demeand;
 subj03_phl_two_demeaned = median_primary_err_hi_left_cond(3,2)- subj03_demeand;
 subj03_phl_three_demeaned = median_primary_err_hi_left_cond(3,3)- subj03_demeand;
 subj03_phr_one_demeaned = median_primary_err_hi_right_cond(3,1)- subj03_demeand;
 subj03_phr_two_demeaned = median_primary_err_hi_right_cond(3,2)- subj03_demeand;
 subj03_phr_three_demeaned = median_primary_err_hi_right_cond(3,3)- subj03_demeand;

 subj03_pll_one_demeaned = median_primary_err_lo_left_cond(3,1)- subj03_demeand;
 subj03_pll_two_demeaned = median_primary_err_lo_left_cond(3,2)- subj03_demeand;
 subj03_pll_three_demeaned = median_primary_err_lo_left_cond(3,3)- subj03_demeand;
 subj03_plr_one_demeaned = median_primary_err_lo_right_cond(3,1)- subj03_demeand;
 subj03_plr_two_demeaned = median_primary_err_lo_right_cond(3,2)- subj03_demeand;
 subj03_plr_three_demeaned = median_primary_err_lo_right_cond(3,3)- subj03_demeand;
  
 subj04_phl_one_demeaned = median_primary_err_hi_left_cond(4,1)- subj04_demeand;
 subj04_phl_two_demeaned = median_primary_err_hi_left_cond(4,2)- subj04_demeand;
 subj04_phl_three_demeaned = median_primary_err_hi_left_cond(4,3)- subj04_demeand;
 subj04_phr_one_demeaned = median_primary_err_hi_right_cond(4,1)- subj04_demeand;
 subj04_phr_two_demeaned = median_primary_err_hi_right_cond(4,2)- subj04_demeand;
 subj04_phr_three_demeaned = median_primary_err_hi_right_cond(4,3)- subj04_demeand;

 subj04_pll_one_demeaned = median_primary_err_lo_left_cond(4,1)- subj04_demeand;
 subj04_pll_two_demeaned = median_primary_err_lo_left_cond(4,2)- subj04_demeand;
 subj04_pll_three_demeaned = median_primary_err_lo_left_cond(4,3)- subj04_demeand;
 subj04_plr_one_demeaned = median_primary_err_lo_right_cond(4,1)- subj04_demeand;
 subj04_plr_two_demeaned = median_primary_err_lo_right_cond(4,2)- subj04_demeand;
 subj04_plr_three_demeaned = median_primary_err_hi_right_cond(4,3)- subj04_demeand;

 
 demeand_pri_hi_left_cond_one = [subj01_phl_one_demeaned subj02_phl_one_demeaned  subj03_phl_one_demeaned subj04_phl_one_demeaned];
 demeand_pri_hi_right_cond_one = [subj01_phr_one_demeaned subj02_phr_one_demeaned  subj03_phr_one_demeaned subj04_phr_one_demeaned];
 demeand_pri_lo_left_cond_one = [subj01_pll_one_demeaned subj02_pll_one_demeaned  subj03_pll_one_demeaned subj04_pll_one_demeaned];
 demeand_pri_lo_right_cond_one = [subj01_plr_one_demeaned subj02_plr_one_demeaned  subj03_plr_one_demeaned subj04_plr_one_demeaned];
 
 % concat group for condition two (l_spcs)
 demeand_pri_hi_left_cond_two = [subj01_phl_two_demeaned subj02_phl_two_demeaned  subj03_phl_two_demeaned subj04_phl_two_demeaned];
 demeand_pri_hi_right_cond_two = [subj01_phr_two_demeaned subj02_phr_two_demeaned  subj03_phr_two_demeaned subj04_phr_two_demeaned];
 demeand_pri_lo_left_cond_two = [subj01_pll_two_demeaned subj02_pll_two_demeaned  subj03_pll_two_demeaned subj04_pll_two_demeaned];
 demeand_pri_lo_right_cond_two = [subj01_plr_two_demeaned subj02_plr_two_demeaned  subj03_plr_two_demeaned subj04_plr_two_demeaned];
 
 demeand_pri_hi_left_cond_three = [subj01_phl_three_demeaned subj02_phl_three_demeaned  subj03_phl_three_demeaned subj04_phl_three_demeaned];
 demeand_pri_hi_right_cond_three = [subj01_phr_three_demeaned subj02_phr_three_demeaned  subj03_phr_three_demeaned subj04_phr_three_demeaned];
 demeand_pri_lo_left_cond_three = [subj01_pll_three_demeaned subj02_pll_three_demeaned  subj03_pll_three_demeaned subj04_pll_three_demeaned];
 demeand_pri_lo_right_cond_three = [subj01_plr_three_demeaned subj02_plr_three_demeaned  subj03_plr_three_demeaned subj04_plr_three_demeaned];
 
ipsi_one =  demeand_pri_lo_left_cond_one - demeand_pri_hi_left_cond_one;
contra_one = demeand_pri_lo_right_cond_one -  demeand_pri_hi_right_cond_one;
ipsi_two =  demeand_pri_lo_left_cond_two -  demeand_pri_hi_left_cond_two;
contra_two = demeand_pri_lo_right_cond_two - demeand_pri_hi_right_cond_two;
ipsi_three =  demeand_pri_lo_left_cond_three - demeand_pri_hi_left_cond_three;
contra_three = demeand_pri_lo_right_cond_three -  demeand_pri_hi_right_cond_three;

sem_ipsi_one = std(ipsi_one)/sqrt(length(ipsi_one));
sem_contra_one = std(contra_one)/sqrt(length(contra_one));
sem_ipsi_two = std(ipsi_two)/sqrt(length(ipsi_two));
sem_contra_two = std(contra_two)/sqrt(length(contra_two));
sem_ipsi_three  = std(ipsi_three )/sqrt(length(ipsi_three));
sem_contra_three  = std(contra_three )/sqrt(length(contra_three ));

demeand_pri_hi_left_cond_one = [subj01_phl_one_demeaned subj02_phl_one_demeaned  subj03_phl_one_demeaned subj04_phl_one_demeaned];
demeand_pri_hi_right_cond_one = [subj01_phr_one_demeaned subj02_phr_one_demeaned  subj03_phr_one_demeaned subj04_phr_one_demeaned];
demeand_pri_lo_left_cond_one = [subj01_pll_one_demeaned subj02_pll_one_demeaned  subj03_pll_one_demeaned subj04_pll_one_demeaned];
demeand_pri_lo_right_cond_one = [subj01_plr_one_demeaned subj02_plr_one_demeaned  subj03_plr_one_demeaned subj04_plr_one_demeaned];
 
 % concat group for condition two - lspcs  
demeand_pri_hi_left_cond_two = [subj01_phl_two_demeaned subj02_phl_two_demeaned  subj03_phl_two_demeaned subj04_phl_two_demeaned];
demeand_pri_hi_right_cond_two = [subj01_phr_two_demeaned subj02_phr_two_demeaned  subj03_phr_two_demeaned subj04_phr_two_demeaned];
demeand_pri_lo_left_cond_two = [subj01_pll_two_demeaned subj02_pll_two_demeaned  subj03_pll_two_demeaned subj04_pll_two_demeaned];
demeand_pri_lo_right_cond_two = [subj01_plr_two_demeaned subj02_plr_two_demeaned  subj03_plr_two_demeaned subj04_plr_two_demeaned];
 
 % concat group for condition three - ips2
demeand_pri_hi_left_cond_three = [subj01_phl_three_demeaned subj02_phl_three_demeaned  subj03_phl_three_demeaned subj04_phl_one_demeaned];
demeand_pri_hi_right_cond_three = [subj01_phr_three_demeaned subj02_phr_three_demeaned  subj03_phr_three_demeaned  subj04_phr_one_demeaned];
demeand_pri_lo_left_cond_three = [subj01_pll_three_demeaned subj02_pll_three_demeaned  subj03_pll_three_demeaned subj04_pll_one_demeaned];
demeand_pri_lo_right_cond_three = [subj01_plr_three_demeaned subj02_plr_three_demeaned  subj03_plr_three_demeaned subj04_plr_one_demeaned];

demeaned_err = [demeand_pri_hi_left_cond_one   demeand_pri_lo_left_cond_one  demeand_pri_hi_right_cond_one demeand_pri_lo_right_cond_one...
demeand_pri_hi_left_cond_two demeand_pri_lo_left_cond_two demeand_pri_hi_right_cond_two demeand_pri_lo_right_cond_two...
demeand_pri_hi_left_cond_three demeand_pri_lo_left_cond_three demeand_pri_hi_right_cond_three demeand_pri_lo_right_cond_three];

diff_vect = [sem_ipsi_one sem_contra_one sem_ipsi_two sem_contra_two sem_ipsi_three sem_contra_three];
%%

hileftsemone = std(demeand_pri_hi_left_cond_one)/sqrt(length(demeand_pri_hi_left_cond_one))
hirightsemone = std(demeand_pri_hi_right_cond_one)/sqrt(length(demeand_pri_hi_right_cond_one))
loleftsemone = std(demeand_pri_lo_left_cond_one)/sqrt(length(demeand_pri_lo_left_cond_one))
lorightsemone = std(demeand_pri_lo_right_cond_one)/sqrt(length(demeand_pri_lo_right_cond_one))

hileftsemtwo = std(demeand_pri_hi_left_cond_two)/sqrt(length(demeand_pri_hi_left_cond_two))
hirightsemtwo= std(demeand_pri_hi_right_cond_two)/sqrt(length(demeand_pri_hi_right_cond_two))
loleftsemtwo = std(demeand_pri_lo_left_cond_two)/sqrt(length(demeand_pri_lo_left_cond_two))
lorightsemtwo = std(demeand_pri_lo_right_cond_two)/sqrt(length(demeand_pri_lo_right_cond_two))

hileftsemthree = std(demeand_pri_hi_left_cond_three)/sqrt(length(demeand_pri_hi_left_cond_three))
hirightsemthree= std(demeand_pri_hi_right_cond_three)/sqrt(length(demeand_pri_hi_right_cond_three))
loleftsemthree = std(demeand_pri_lo_left_cond_three)/sqrt(length(demeand_pri_lo_left_cond_three))
lorightsemthree = std(demeand_pri_lo_right_cond_three)/sqrt(length(demeand_pri_lo_right_cond_three))

 sem_vect= [hileftsemone  loleftsemone hirightsemone lorightsemone...
 hileftsemtwo loleftsemtwo hirightsemtwo lorightsemtwo...
 hileftsemthree loleftsemthree hirightsemthree lorightsemthree];

%% plotting loop - PRIMARY, demeaned

 for ii=1:3
figure(ii)
subplot(1,3,1)
hold on;
errorbar([demeand_pri_hi_left_cond_one(ii) demeand_pri_lo_left_cond_one(ii)] , [sem_vect(1)  sem_vect(2)] ,'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1)%'.')
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
set(gca, 'fontsize', 14)
hold on; 
errorbar([demeand_pri_hi_right_cond_one(ii) demeand_pri_lo_right_cond_one(ii)], [sem_vect(3) sem_vect(4)] ,'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1)%'.')
legend('No TMS ipsi', 'No TMS contra')

subplot(1,3,2)
hold on;
errorbar([demeand_pri_hi_left_cond_two(ii) demeand_pri_lo_left_cond_two(ii)] , [sem_vect(5) sem_vect(6)] ,'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13)
errorbar([demeand_pri_hi_right_cond_two(ii) demeand_pri_lo_right_cond_two(ii)] , [sem_vect(7) sem_vect(8)] ,'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13)
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
set(gca, 'fontsize', 14)
legend('sPCS ipsi','sPCS contra')
hold on; 

subplot(1,3,3)
errorbar([demeand_pri_hi_left_cond_three(ii) demeand_pri_lo_left_cond_three(ii)] , [sem_vect(9) sem_vect(10)] ,'o-','color', c5,'markersize',10,'MarkerEdgeColor', c5)
hold on;
errorbar([demeand_pri_hi_right_cond_three(ii) demeand_pri_lo_right_cond_three(ii)] , [sem_vect(11) sem_vect(12)] ,'o--','color', c5,'markersize', 10,'MarkerEdgeColor', c5,'MarkerFaceColor', c5)
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
legend('IPS2 ipsi','IPS2 contra')
set(gca, 'fontsize', 14)
 end 
%% PRIMARY - demeaned error bars only
 for ii =1:2;
group_mean_lo_right(ii) = mean(median_primary_err_lo_right_cond(:,ii)); %takes mean down columns (aka, each row is a subj, so mean of all subj in that cond)
group_mean_lo_left(ii) = mean(median_primary_err_lo_left_cond(:,ii));
group_mean_hi_left(ii) = mean(median_primary_err_hi_left_cond(:,ii));
group_mean_hi_right(ii) = mean(median_primary_err_hi_right_cond(:,ii));
end 

for ii = 1;
    figure(ii)
    hold on;
    errorbar([group_mean_hi_left(ii) group_mean_lo_left(ii)] , [sem_vect(1)  sem_vect(2)] ,'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1)%'.')
    set(gca, 'xtick', [1 2])
    set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
    set(gca, 'fontsize', 14)
    ylim([0.5 2.5])
    hold on;
    errorbar([group_mean_hi_right(ii) group_mean_lo_right(ii)], [sem_vect(3) sem_vect(4)] ,'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1)%'.')
    legend('No TMS ipsi', 'No TMS contra')
end

for jj = 2
    figure(jj)
    errorbar([group_mean_hi_left(jj) group_mean_lo_left(jj)] , [sem_vect(5) sem_vect(6)] ,'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13)
    hold on;
    errorbar([group_mean_hi_right(jj) group_mean_lo_right(jj)] , [sem_vect(7) sem_vect(8)] ,'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13)
    set(gca, 'xtick', [1 2])
    set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
    set(gca, 'fontsize', 14)
    ylim([0.5 2.5])
    legend('sPCS ipsi','sPCS contra')
end

 for kk = 3
figure(kk)
errorbar([group_mean_hi_left(kk) group_mean_lo_left(kk)] , [sem_vect(9) sem_vect(10)] ,'o-','color', c5,'markersize',10,'MarkerEdgeColor', c5)
hold on;
errorbar([group_mean_hi_right(kk) group_mean_lo_right(kk)] , [sem_vect(11)  sem_vect(12)] ,'o--','color', c5,'markersize', 10,'MarkerEdgeColor', c5,'MarkerFaceColor', c5)
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
legend('IPS2 ipsi','IPS2 contra')
set(gca, 'fontsize', 14)
ylim([0.5 2.5])
 end
 
 %% PRIMARY, subject-wise avg with demeaned error bars
% subj01_demeand = mean([subj_mean_pri_hi_left(1) subj_mean_pri_hi_right(1) subj_mean_pri_lo_left(1) subj_mean_pri_lo_right(1)]);
% subj02_demeand = mean([subj_mean_pri_hi_left(2) subj_mean_pri_hi_right(2) subj_mean_pri_lo_left(2) subj_mean_pri_lo_right(2)]); 
% subj03_demeand = mean([subj_mean_pri_hi_left(3) subj_mean_pri_hi_right(3) subj_mean_pri_lo_left(3) subj_mean_pri_lo_right(3)]);
% subj04_demeand = mean([subj_mean_pri_hi_left(4) subj_mean_pri_hi_right(4) subj_mean_pri_lo_left(4) subj_mean_pri_lo_right(4)]);


 for ii=1
figure(ii)
subplot(1,3,1)
hold on;
errorbar([subj_mean_pri_hi_left(ii) subj_mean_pri_lo_left(ii)] , [sem_vect(1)  sem_vect(2)] ,'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1)%'.')
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
set(gca, 'fontsize', 14)
hold on; 
errorbar([subj_mean_pri_hi_right(ii) subj_mean_pri_lo_right(ii)], [sem_vect(3) sem_vect(4)] ,'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1)%'.')
legend('No TMS ipsi', 'No TMS contra')

subplot(1,3,2)
hold on;
errorbar([subj_mean_pri_hi_left(ii) subj_mean_pri_lo_left(ii)] , [sem_vect(5) sem_vect(6)] ,'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13)
errorbar([subj_mean_pri_hi_right(ii) subj_mean_pri_lo_right(ii)] , [sem_vect(7) sem_vect(8)] ,'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13)
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
set(gca, 'fontsize', 14)
legend('sPCS ipsi','sPCS contra')
hold on; 

subplot(1,3,3)
errorbar([subj_mean_pri_hi_left(ii) subj_mean_pri_lo_left(ii)] , [sem_vect(9) sem_vect(10)] ,'o-','color', c5,'markersize',10,'MarkerEdgeColor', c5)
hold on;
errorbar([subj_mean_pri_hi_right(ii) subj_mean_pri_lo_right(ii)] , [sem_vect(11) sem_vect(12)] ,'o--','color', c5,'markersize', 10,'MarkerEdgeColor', c5,'MarkerFaceColor', c5)
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
legend('IPS2 ipsi','IPS2 contra')
set(gca, 'fontsize', 14)
 end 
 
 
 
 
 
 
 %%  FINAL

subj = {'subj01','subj02','subj03','subj04'};
cond  =  {'pilot','l_ips2','l_spcs'}; %'l_ips2',
num_cond = length(cond);
num_subj = length(subj);  


for ss = 1:num_subj;
    for cc = 1:num_cond;
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
        resultsfile = load(filename);
        median_final_err_lo_left_cond(ss,cc) = resultsfile.ii_results_lo.median_no_break_left_final_err_z_new;
        median_final_err_lo_right_cond(ss,cc) = resultsfile.ii_results_lo.median_no_break_right_final_err_z_new;
    end 
    for cc = 1:num_cond;
        filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
        resultsfile = load(filename)
        median_final_err_hi_left_cond(ss,cc) = resultsfile.ii_results_hi.median_no_break_left_final_err_z_new;
        median_final_err_hi_right_cond(ss,cc) = resultsfile.ii_results_hi.median_no_break_right_final_err_z_new;
    end
end 

for ii = 1:num_cond;
    group_mean_hi_left(ii) = mean(median_final_err_hi_left_cond(:,ii));
    group_mean_hi_right(ii) = mean(median_final_err_hi_right_cond(:,ii));
    group_mean_lo_right(ii) = mean(median_final_err_lo_right_cond(:,ii));
    group_mean_lo_left(ii) = mean(median_final_err_lo_left_cond(:,ii));
end

for ii = 1:num_cond;
    all_cond = [];
    all_cond = [group_mean_hi_left group_mean_hi_right group_mean_lo_left group_mean_lo_right;];
    
end

for jj = 1:num_subj;
    subj_mean_final_hi_left = mean(median_final_err_hi_left_cond,2)
    subj_mean_final_hi_right = mean(median_final_err_hi_right_cond,2)
    subj_mean_final_lo_left = mean(median_final_err_lo_left_cond,2)
    subj_mean_final_lo_right = mean(median_final_err_lo_right_cond,2)
end

%% 
subj01_demeand = mean([subj_mean_final_hi_left(1) subj_mean_final_hi_right(1) subj_mean_final_lo_left(1) subj_mean_final_lo_right(1)]);
subj02_demeand = mean([subj_mean_final_hi_left(2) subj_mean_final_hi_right(2) subj_mean_final_lo_left(2) subj_mean_final_lo_right(2)]); 
subj03_demeand = mean([subj_mean_final_hi_left(3) subj_mean_final_hi_right(3) subj_mean_final_lo_left(3) subj_mean_final_lo_right(3)]);
subj04_demeand = mean([subj_mean_final_hi_left(4) subj_mean_final_hi_right(4) subj_mean_final_lo_left(4) subj_mean_final_lo_right(4)]);


 subj01_phl_one_demeaned = median_final_err_hi_left_cond(1,1)- subj01_demeand;
 subj01_phl_two_demeaned = median_final_err_hi_left_cond(1,2)- subj01_demeand;
 subj01_phl_three_demeaned = median_final_err_hi_left_cond(1,3)- subj01_demeand;
 subj01_phr_one_demeaned = median_final_err_hi_right_cond(1,1)- subj01_demeand;
 subj01_phr_two_demeaned = median_final_err_hi_right_cond(1,2)- subj01_demeand;
 subj01_phr_three_demeaned = median_final_err_hi_right_cond(1,3)- subj01_demeand;

 subj01_pll_one_demeaned = median_final_err_lo_left_cond(1,1)- subj01_demeand;
 subj01_pll_two_demeaned = median_final_err_lo_left_cond(1,2)- subj01_demeand;
 subj01_pll_three_demeaned = median_final_err_lo_left_cond(1,3)- subj01_demeand;
 subj01_plr_one_demeaned = median_final_err_lo_right_cond(1,1)- subj01_demeand;
 subj01_plr_two_demeaned = median_final_err_lo_right_cond(1,2)- subj01_demeand;
 subj01_plr_three_demeaned = median_final_err_lo_right_cond(1,3)- subj01_demeand;
 
 
 subj02_phl_one_demeaned = median_final_err_hi_left_cond(2,1)- subj02_demeand;
 subj02_phl_two_demeaned = median_final_err_hi_left_cond(2,2)- subj02_demeand;
 subj02_phl_three_demeaned = median_final_err_hi_left_cond(2,3)- subj02_demeand;
 subj02_phr_one_demeaned = median_final_err_hi_right_cond(2,1)- subj02_demeand;
 subj02_phr_two_demeaned = median_final_err_hi_right_cond(2,2)- subj02_demeand;
 subj02_phr_three_demeaned = median_final_err_hi_right_cond(2,3)- subj02_demeand;

 subj02_pll_one_demeaned = median_final_err_lo_left_cond(2,1)- subj02_demeand;
 subj02_pll_two_demeaned = median_final_err_lo_left_cond(2,2)- subj02_demeand;
 subj02_pll_three_demeaned = median_final_err_lo_left_cond(2,3)- subj02_demeand;
 subj02_plr_one_demeaned = median_final_err_lo_right_cond(2,1)- subj02_demeand;
 subj02_plr_two_demeaned = median_final_err_lo_right_cond(2,2)- subj02_demeand;
 subj02_plr_three_demeaned = median_final_err_lo_right_cond(2,3)- subj02_demeand;
 
 subj03_phl_one_demeaned = median_final_err_hi_left_cond(3,1)- subj03_demeand;
 subj03_phl_two_demeaned = median_final_err_hi_left_cond(3,2)- subj03_demeand;
 subj03_phl_three_demeaned = median_final_err_hi_left_cond(3,3)- subj03_demeand;
 subj03_phr_one_demeaned = median_final_err_hi_right_cond(3,1)- subj03_demeand;
 subj03_phr_two_demeaned = median_final_err_hi_right_cond(3,2)- subj03_demeand;
 subj03_phr_three_demeaned = median_final_err_hi_right_cond(3,3)- subj03_demeand;

 subj03_pll_one_demeaned = median_final_err_lo_left_cond(3,1)- subj03_demeand;
 subj03_pll_two_demeaned = median_final_err_lo_left_cond(3,2)- subj03_demeand;
 subj03_pll_three_demeaned = median_final_err_lo_left_cond(3,3)- subj03_demeand;
 subj03_plr_one_demeaned = median_final_err_lo_right_cond(3,1)- subj03_demeand;
 subj03_plr_two_demeaned = median_final_err_lo_right_cond(3,2)- subj03_demeand;
 subj03_plr_three_demeaned = median_final_err_lo_right_cond(3,3)- subj03_demeand;
 %  
 subj04_phl_one_demeaned = median_final_err_hi_left_cond(4,1)- subj04_demeand;
 subj04_phl_two_demeaned = median_final_err_hi_left_cond(4,2)- subj04_demeand;
 subj04_phl_three_demeaned = median_final_err_hi_left_cond(4,3)- subj04_demeand;
 subj04_phr_one_demeaned = median_final_err_hi_right_cond(4,1)- subj04_demeand;
 subj04_phr_two_demeaned = median_final_err_hi_right_cond(4,2)- subj04_demeand;
 subj04_phr_three_demeaned = median_final_err_hi_right_cond(4,3)- subj04_demeand;

 subj04_pll_one_demeaned = median_final_err_lo_left_cond(4,1)- subj04_demeand;
 subj04_pll_two_demeaned = median_final_err_lo_left_cond(4,2)- subj04_demeand;
 subj04_pll_three_demeaned = median_final_err_lo_left_cond(4,3)- subj04_demeand;
 subj04_plr_one_demeaned = median_final_err_lo_right_cond(4,1)- subj04_demeand;
 subj04_plr_two_demeaned = median_final_err_lo_right_cond(4,2)- subj04_demeand;
 subj04_plr_three_demeaned = median_final_err_hi_right_cond(4,3)- subj04_demeand;

 
 demeand_pri_hi_left_cond_one = [subj01_phl_one_demeaned subj02_phl_one_demeaned  subj03_phl_one_demeaned subj04_phl_one_demeaned];
 demeand_pri_hi_right_cond_one = [subj01_phr_one_demeaned subj02_phr_one_demeaned  subj03_phr_one_demeaned subj04_phr_one_demeaned];
 demeand_pri_lo_left_cond_one = [subj01_pll_one_demeaned subj02_pll_one_demeaned  subj03_pll_one_demeaned  subj04_pll_one_demeaned];
 demeand_pri_lo_right_cond_one = [subj01_plr_one_demeaned subj02_plr_one_demeaned  subj03_plr_one_demeaned  subj04_plr_one_demeaned];
 
 % concat group for condition two (l_spcs)
 demeand_pri_hi_left_cond_two = [subj01_phl_two_demeaned subj02_phl_two_demeaned  subj03_phl_two_demeaned subj04_phl_two_demeaned];
 demeand_pri_hi_right_cond_two = [subj01_phr_two_demeaned subj02_phr_two_demeaned  subj03_phr_two_demeaned subj04_phr_two_demeaned];
 demeand_pri_lo_left_cond_two = [subj01_pll_two_demeaned subj02_pll_two_demeaned  subj03_pll_two_demeaned subj04_pll_two_demeaned];
 demeand_pri_lo_right_cond_two = [subj01_plr_two_demeaned subj02_plr_two_demeaned  subj03_plr_two_demeaned subj04_plr_two_demeaned];
 
 demeand_pri_hi_left_cond_three = [subj01_phl_three_demeaned subj02_phl_three_demeaned  subj03_phl_three_demeaned subj04_phl_three_demeaned];
 demeand_pri_hi_right_cond_three = [subj01_phr_three_demeaned subj02_phr_three_demeaned  subj03_phr_three_demeaned subj04_phr_three_demeaned];
 demeand_pri_lo_left_cond_three = [subj01_pll_three_demeaned subj02_pll_three_demeaned  subj03_pll_three_demeaned subj04_pll_three_demeaned];
 demeand_pri_lo_right_cond_three = [subj01_plr_three_demeaned subj02_plr_three_demeaned  subj03_plr_three_demeaned subj04_plr_three_demeaned];
 
ipsi_one =  demeand_pri_lo_left_cond_one - demeand_pri_hi_left_cond_one;
contra_one = demeand_pri_lo_right_cond_one -  demeand_pri_hi_right_cond_one;
ipsi_two =  demeand_pri_lo_left_cond_two -  demeand_pri_hi_left_cond_two;
contra_two = demeand_pri_lo_right_cond_two - demeand_pri_hi_right_cond_two;
ipsi_three =  demeand_pri_lo_left_cond_three - demeand_pri_hi_left_cond_three;
contra_three = demeand_pri_lo_right_cond_three -  demeand_pri_hi_right_cond_three;

sem_ipsi_one = std(ipsi_one)/sqrt(length(ipsi_one));
sem_contra_one = std(contra_one)/sqrt(length(contra_one));
sem_ipsi_two = std(ipsi_two)/sqrt(length(ipsi_two));
sem_contra_two = std(contra_two)/sqrt(length(contra_two));
sem_ipsi_three  = std(ipsi_three )/sqrt(length(ipsi_three));
sem_contra_three  = std(contra_three )/sqrt(length(contra_three ));

demeand_final_hi_left_cond_one = [subj01_phl_one_demeaned subj02_phl_one_demeaned  subj03_phl_one_demeaned subj04_phl_one_demeaned];
demeand_final_hi_right_cond_one = [subj01_phr_one_demeaned subj02_phr_one_demeaned  subj03_phr_one_demeaned subj04_phr_one_demeaned];
demeand_final_lo_left_cond_one = [subj01_pll_one_demeaned subj02_pll_one_demeaned  subj03_pll_one_demeaned subj04_pll_one_demeaned];
demeand_final_lo_right_cond_one = [subj01_plr_one_demeaned subj02_plr_one_demeaned  subj03_plr_one_demeaned subj04_plr_one_demeaned];
 
 % concat group for condition two - lspcs  
demeand_final_hi_left_cond_two = [subj01_phl_two_demeaned subj02_phl_two_demeaned  subj03_phl_two_demeaned subj04_phl_two_demeaned];
demeand_final_hi_right_cond_two = [subj01_phr_two_demeaned subj02_phr_two_demeaned  subj03_phr_two_demeaned subj04_phr_two_demeaned];
demeand_final_lo_left_cond_two = [subj01_pll_two_demeaned subj02_pll_two_demeaned  subj03_pll_two_demeaned subj04_pll_two_demeaned];
demeand_final_lo_right_cond_two = [subj01_plr_two_demeaned subj02_plr_two_demeaned  subj03_plr_two_demeaned subj04_plr_two_demeaned];
 
 % concat group for condition three - ips2
demeand_final_hi_left_cond_three = [subj01_phl_three_demeaned subj02_phl_three_demeaned  subj03_phl_three_demeaned subj04_phl_one_demeaned];
demeand_final_hi_right_cond_three = [subj01_phr_three_demeaned subj02_phr_three_demeaned  subj03_phr_three_demeaned subj04_phr_one_demeaned];
demeand_final_lo_left_cond_three = [subj01_pll_three_demeaned subj02_pll_three_demeaned  subj03_pll_three_demeaned subj04_pll_one_demeaned];
demeand_final_lo_right_cond_three = [subj01_plr_three_demeaned subj02_plr_three_demeaned  subj03_plr_three_demeaned subj04_plr_one_demeaned];

demeaned_err = [demeand_final_hi_left_cond_one demeand_final_lo_left_cond_one  demeand_final_hi_right_cond_one demeand_final_lo_right_cond_one...
demeand_final_hi_left_cond_two demeand_pri_lo_left_cond_two demeand_final_hi_right_cond_two demeand_final_lo_right_cond_two...
demeand_final_hi_left_cond_three demeand_final_lo_left_cond_three demeand_final_hi_right_cond_three demeand_final_lo_right_cond_three];

diff_vect = [sem_ipsi_one sem_contra_one sem_ipsi_two sem_contra_two sem_ipsi_three sem_contra_three];

%%  for ii= 1:4

for ii = 1:4
figure(ii);
subplot(1,3,1)
hold on;
plot([1,2], [demeand_final_hi_left_cond_one(ii)  demeand_final_lo_left_cond_one(ii)],'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1) %ipsi
errorbar([demeand_final_hi_left_cond_one(ii) demeand_final_lo_left_cond_one(ii)] , [diff_vect(1) diff_vect(1)] ,'.')
plot([1,2], [demeand_final_hi_right_cond_one(ii)  demeand_final_lo_right_cond_one(ii)],'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1) %contra
hold on; 
errorbar([demeand_final_hi_right_cond_one(ii) demeand_final_lo_right_cond_one(ii)], [diff_vect(2) diff_vect(2)] ,'.')
ylim([-1 1])
subplot(1,3,2)
plot([1,2], [demeand_final_hi_left_cond_two(ii)  demeand_final_lo_left_cond_two(ii)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13) %ipsi
hold on;
errorbar([demeand_final_hi_left_cond_two(ii) demeand_final_lo_left_cond_two(ii)] , [diff_vect(3) diff_vect(3)] ,'.')
plot([1,2], [demeand_final_hi_right_cond_two(ii)  demeand_final_lo_right_cond_two(ii)],'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13) %contra
errorbar([demeand_final_hi_right_cond_two(ii) demeand_final_lo_right_cond_two(ii)] , [diff_vect(4) diff_vect(4)] ,'.')
ylim([-1 1])
subplot(1,3,3)
plot([1,2], [demeand_final_hi_left_cond_three(ii)  demeand_final_lo_left_cond_three(ii)],'o-','color', c5,'markersize',10,'MarkerEdgeColor', c5) %ipsi
hold on;
errorbar([demeand_final_hi_left_cond_three(ii) demeand_final_lo_left_cond_three(ii)] , [diff_vect(5) diff_vect(5)] ,'.')
plot([1,2], [demeand_final_hi_right_cond_three(ii)  demeand_final_lo_right_cond_three(ii)],'o--','color', c5,'markersize', 10,'MarkerEdgeColor', c5,'MarkerFaceColor', c5) %contra
errorbar([demeand_final_hi_right_cond_three(ii) demeand_final_lo_right_cond_three(ii)] , [diff_vect(6) diff_vect(6)] ,'.')
ylim([-1 1])
 end 
