%% make plots of primary and final eye positions
%% custom colors
par = parula;
beta = 0.1;
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

% c11 = par(39,:,:)
% c12 = par(40,:,:)

c13 = par(55,:,:);
c14 = par(56,:,:);
c15 = brighten(c13,beta);
c16 = brighten(c14,beta);



% c15 = par(63,:,:)
% c16 = par(64,:,:)
mycolors_i = {c1 c13 c9};
mycolors_c = {c3 c15 c11};
parspec = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16];
%% PRIMARY EYE POSITION
%subject = {'subj01', 'subj02', 'subj03','subj04'};
subject = {'subj07'};

%subj = {'gh','JF','MP','pk','EK','mr','cc'};
%subj = {'gh','MP','pk','EK','mr','cc'};

for jj =1:length(subject);
subj = {[sprintf('%s',subject{jj})']};
%cond = {'noTMS','l_spcs','l_ips2'};
cond={'noTMS'};
primary_err_lo_left_subj = [];
median_primary_err_lo_left_group = [];
primary_err_lo_left_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    primary_err_lo_left_subj{ss,cc} =  resultsfile.ii_results_lo.no_break_left_primary_err_z_new;
    median_primary_err_lo_left(ss,cc) = [median_primary_err_lo_left_group;
    resultsfile.ii_results_lo.median_no_break_left_primary_err_z_new];
    primary_err_lo_left_sem(ss,cc) = std(primary_err_lo_left_subj{ss,cc})/sqrt(length(primary_err_lo_left_subj{ss,cc}))
 
end
end

subj = {[sprintf('%s',subject{jj})']};
%cond = {'noTMS','l_spcs','l_ips2'};
cond={'noTMS'};
primary_err_lo_right_subj = [];
median_primary_err_lo_right_group = [];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    primary_err_lo_right_subj{ss,cc} =  resultsfile.ii_results_lo.no_break_right_primary_err_z_new;
    median_primary_err_lo_right(ss,cc) = [median_primary_err_lo_right_group;
    resultsfile.ii_results_lo.median_no_break_right_primary_err_z_new];
    primary_err_lo_right_sem(ss,cc) = std(primary_err_lo_right_subj{ss,cc})/sqrt(length(primary_err_lo_right_subj{ss,cc}))
 
end
end


subj = {[sprintf('%s',subject{jj})']};
%cond = {'noTMS','l_spcs','l_ips2'};
cond={'noTMS'};
primary_err_hi_left_subj = [];
median_primary_err_hi_left_group = [];
primary_err_hi_left_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    primary_err_hi_left_subj{ss,cc} =  resultsfile.ii_results_hi.no_break_left_primary_err_z_new;
    median_primary_err_hi_left(ss,cc) = [median_primary_err_hi_left_group;
    resultsfile.ii_results_hi.median_no_break_left_primary_err_z_new];
    primary_err_hi_left_sem(ss,cc) = std(primary_err_hi_left_subj{ss,cc})/sqrt(length(primary_err_hi_left_subj{ss,cc}))
 
end
end

subj = {[sprintf('%s',subject{jj})']};
%cond = {'noTMS','l_spcs','l_ips2'};
cond={'noTMS'};
primary_err_hi_right_subj = [];
median_primary_err_hi_right_group = [];
primary_err_hi_right_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    primary_err_hi_right_subj{ss,cc} =  resultsfile.ii_results_hi.no_break_right_primary_err_z_new;
    median_primary_err_hi_right(ss,cc) = [median_primary_err_hi_right_group;
    resultsfile.ii_results_hi.median_no_break_right_primary_err_z_new];
    primary_err_hi_right_sem(ss,cc) = std(primary_err_hi_right_subj{ss,cc})/sqrt(length(primary_err_hi_right_subj{ss,cc}))
 
end
end


%% plot the above 
subj = {[sprintf('%s',subject{jj})']};
%cond = {'noTMS','l_spcs','l_ips2'};
cond={'noTMS'};
for cc =1:length(cond);
    figure(jj);
    errorbar([median_primary_err_hi_left(ss,cc)  median_primary_err_lo_left(ss,cc)], [primary_err_hi_left_sem(ss,cc) primary_err_lo_left_sem(ss,cc)],'o-','color', mycolors_i{cc},'markersize', 10,'MarkerEdgeColor', mycolors_i{cc});
    hold on;
    legend('no tms ipsi', 'spcs ipsi', 'ips ipsi', 'no tms contra', 'spcs contra', 'ips contra')
    [vect1 h] = ranksum(primary_err_hi_left_subj{cc}, primary_err_lo_left_subj{cc})
    if vect1 < 0.05
        h1  = sigstar([1 2],vect1)
        set(h1,'color',mycolors_i{cc})
    else
    end
    errorbar([median_primary_err_hi_right(ss,cc)  median_primary_err_lo_right(ss,cc)], [primary_err_hi_right_sem(ss,cc) primary_err_lo_right_sem(ss,cc)],'o--','color', mycolors_c{cc},'markersize', 10,'MarkerEdgeColor', mycolors_c{cc},'MarkerFaceColor', mycolors_c{cc});
    hold on;
      %legend([sprintf('%s',cond{cc}),' contra'])
    [vect3 h] = ranksum(primary_err_hi_right_subj{cc}, primary_err_lo_right_subj{cc})
    if vect3 < 0.05
        h2  = sigstar([1 2],vect3)
        set(h2,'color',mycolors_c{cc})
        set(h2,'LineStyle','--')
    else
    end
     
end
ylabel('Absolute Error, DVA')
xlabel('Priority')
title(['Primary',sprintf(' %s',subj{ss})])
%legend('No TMS ipsi','No TMS contra', 'sPCS ipsi','sPCS contra','IPS2 ipsi','IPS2 contra')
set(gca,'FontSize',14)
set(gca,'Xtick', [1 2])
ylim([0.5 3])
set(gca,'Xticklabel', {'High','Low'})

end



[vectspcscontra h] = ranksum(primary_err_hi_right_subj{1}, primary_err_hi_right_subj{2})


%%
% 
% [vect1 h] = ranksum(primary_err_hi_left_subj, primary_err_lo_left_subj)
% [vect2 h] = ranksum(primary_err_hi_left_subj, primary_err_hi_right_subj)
% [vect3 h] = ranksum(primary_err_hi_right_subj, primary_err_lo_right_subj)
% [vect4 h] = ranksum(primary_err_lo_left_subj, primary_err_lo_right_subj)
% 
% sigvect = [vect1 vect3];
% sigstar({[1,3],[2,4]},sigvect)


%% final EYE POSITION
subject = {'subj01', 'subj02', 'subj03','subj04'};
%subject = {'subj04'};
%need to be able to access means of each subj 
%subj = {'gh','JF','MP','pk','EK','mr','cc'};
%subj = {'gh','MP','pk','EK','mr','cc'};

for jj =1:length(subject);
subj = {[sprintf('%s',subject{jj})']};
cond = {'noTMS','l_spcs','l_ips2'};
final_err_lo_left_subj = [];
median_final_err_lo_left_group = [];
final_err_lo_left_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    final_err_lo_left_subj{ss,cc} =  resultsfile.ii_results_lo.no_break_left_final_err_z_new;
    median_final_err_lo_left(ss,cc) = [median_final_err_lo_left_group;
    resultsfile.ii_results_lo.median_no_break_left_final_err_z_new];
    final_err_lo_left_sem(ss,cc) = std(final_err_lo_left_subj{ss,cc})/sqrt(length(final_err_lo_left_subj{ss,cc}))
 
end
end

subj = {[sprintf('%s',subject{jj})']};
cond = {'noTMS','l_spcs','l_ips2'};
final_err_lo_right_subj = [];
median_final_err_lo_right_group = [];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    final_err_lo_right_subj{ss,cc} =  resultsfile.ii_results_lo.no_break_right_final_err_z_new;
    median_final_err_lo_right(ss,cc) = [median_final_err_lo_right_group;
    resultsfile.ii_results_lo.median_no_break_right_final_err_z_new];
    final_err_lo_right_sem(ss,cc) = std(final_err_lo_right_subj{ss,cc})/sqrt(length(final_err_lo_right_subj{ss,cc}))
 
end
end


subj = {[sprintf('%s',subject{jj})']};
cond = {'noTMS','l_spcs','l_ips2'};
final_err_hi_left_subj = [];
median_final_err_hi_left_group = [];
final_err_hi_left_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    final_err_hi_left_subj{ss,cc} =  resultsfile.ii_results_hi.no_break_left_final_err_z_new;
    median_final_err_hi_left(ss,cc) = [median_final_err_hi_left_group;
    resultsfile.ii_results_hi.median_no_break_left_final_err_z_new];
    final_err_hi_left_sem(ss,cc) = std(final_err_hi_left_subj{ss,cc})/sqrt(length(final_err_hi_left_subj{ss,cc}))
 
end
end

subj = {[sprintf('%s',subject{jj})']};
cond = {'noTMS','l_spcs','l_ips2'};
final_err_hi_right_subj = [];
median_final_err_hi_right_group = [];
final_err_hi_right_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    final_err_hi_right_subj{ss,cc} =  resultsfile.ii_results_hi.no_break_right_final_err_z_new;
    median_final_err_hi_right(ss,cc) = [median_final_err_hi_right_group;
    resultsfile.ii_results_hi.median_no_break_right_final_err_z_new];
    final_err_hi_right_sem(ss,cc) = std(final_err_hi_right_subj{ss,cc})/sqrt(length(final_err_hi_right_subj{ss,cc}))
 
end
end


%% plot the above 
subj = {[sprintf('%s',subject{jj})']};
cond = {'noTMS','l_spcs','l_ips2'};
for cc =1:length(cond);
    figure(jj);
    errorbar([median_final_err_hi_left(ss,cc)  median_final_err_lo_left(ss,cc)], [final_err_hi_left_sem(ss,cc) final_err_lo_left_sem(ss,cc)],'o-','color', mycolors_i{cc},'markersize', 10,'MarkerEdgeColor', mycolors_i{cc});
    hold on;
    legend('no tms ipsi', 'spcs ipsi', 'ips ipsi', 'no tms contra', 'spcs contra', 'ips contra')
    [vect1 h] = ranksum(final_err_hi_left_subj{cc}, final_err_lo_left_subj{cc})
    if vect1 < 0.05
        h1  = sigstar([1 2],vect1)
        set(h1,'color',mycolors_i{cc})
    else
    end
    errorbar([median_final_err_hi_right(ss,cc)  median_final_err_lo_right(ss,cc)], [final_err_hi_right_sem(ss,cc) final_err_lo_right_sem(ss,cc)],'o--','color', mycolors_c{cc},'markersize', 10,'MarkerEdgeColor', mycolors_c{cc},'MarkerFaceColor', mycolors_c{cc});
    hold on;
      %legend([sprintf('%s',cond{cc}),' contra'])
    [vect3 h] = ranksum(final_err_hi_right_subj{cc}, final_err_lo_right_subj{cc})
    if vect3 < 0.05
        h2  = sigstar([1 2],vect3)
        set(h2,'color',mycolors_c{cc})
        set(h2,'LineStyle','--')
    else
    end
     
end
ylabel('Absolute Error, DVA')
xlabel('Priority')
title(['final',sprintf(' %s',subj{ss})])
%legend('No TMS ipsi','No TMS contra', 'sPCS ipsi','sPCS contra','IPS2 ipsi','IPS2 contra')
set(gca,'FontSize',14)
set(gca,'Xtick', [1 2])
ylim([0.5 1.75])
set(gca,'Xticklabel', {'High','Low'})

end


%%  coding experiment
%loop over ii_results.files & store them in a single structure
%goal: all subj embedded in one struct we can refer to with descriptive
%varaible names 

subj = {'gh','JF','MP','pk','EK','mr','cc'};
hemi = {'left','right'}; %left or right hemifield 
alldata = [];
cond = {'hi','lo'}; %priority cond
em = {'primary', 'final'}; %eye movement cat
barvect_hi = [];
barvect_lo = [];

        for ss = 1:length(subj);
            for cc = 1:length(cond);
                for hh = 1:length(hemi);
                    for ee = 1:length(em);
                        datafn = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_%s',subj{ss},cond{cc});
                        %data = load(datafn,sprintf('ii_results_%s',cond{cc}));
                        data = load(datafn,sprintf('ii_results_%s',cond{cc}));
                        %need to be in ii_results_hi here in order to refer
                        %within it in the next line like so: 'ii_results_hi.no_break_left_primary_err_z'
                        alldata.(subj{ss}).(cond{cc}).(hemi{hh}).(em{ee}) = data.(sprintf('ii_results_%s',cond{cc})).(sprintf('no_break_%s_%s_err_z',hemi{hh},em{ee})); 
                      %alldata.(subj{ss}).(cond{cc}).(hemi{hh}).(em{ee}) = data.(sprintf('ii_results_%s.no_break_%s_%s_err_z',cond{cc},hemi{hh},em{ee})); %wont load. 
                       %alldatatest.(subj{ss}).(cond{cc}).(hemi{hh}).(em{ee}) = load('alldata');
                       allmeds.(subj{ss}).(cond{cc}).(hemi{hh}).(em{ee}) = data.ii_results_hi.(sprintf('median_no_break_%s_%s_err_z',hemi{hh},em{ee})); 
                    end
                end
                
            end
        end
 
 % stopping point 6/12-- having trouble loading incorporating both
 % 'ii_results_hi' & 'ii_results_lo' from subj folders. can't sprintf into
 % both in "data" step bc not sure how to call it from data.ii_results_??
 % in the next step. 


for jj =  1:length(subj);
barvect_hi_pri(jj,:) = [allmeds.(subj{jj}).hi.left.primary;  allmeds.(subj{jj}).hi.right.primary];
barvect_hi_final(jj,:) = [allmeds.(subj{jj}).hi.left.final;  allmeds.(subj{jj}).hi.right.final];
figure(jj);
bar([barvect_hi_pri(jj,1) barvect_hi_pri(jj,2)]);
end



    
    



