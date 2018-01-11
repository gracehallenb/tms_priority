%% make plots of primary and final eye positions
%% colors
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
%% PRIMARY EYE POSITION
%need to be able to access means of each subj 
%subj = {'gh','JF','MP','pk','EK','mr','cc'};
%subj = {'gh','MP','pk','EK','mr','cc'};
subj = {'subj02'};
%subj = {'subj01','subj02'};
cond = {'noTMS'};
primary_err_lo_left_subj = [];
median_primary_err_lo_left_group = [];
primary_err_lo_left_subj_ind =[];

for cc =1:length(cond);
for ss = 1:length(subj);
    %filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_lo.mat',subj{ss});
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename)
    primary_err_lo_left_subj =  [primary_err_lo_left_subj;   resultsfile.ii_results_lo.no_break_left_primary_err_z];
    %primary_err_lo_left_subj{ss} =  resultsfile.ii_results_lo.no_break_left_primary_err_z;
    %median_primary_err_lo_left_group = [median_primary_err_lo_left_group;
    %resultsfile.ii_results_lo.median_no_break_left_primary_err_z];
    %%commented out 9/13
     median_primary_err_lo_left_group = [median_primary_err_lo_left_group; resultsfile.ii_results_lo.median_no_break_left_primary_err_z];
    primary_err_lo_left_group_mean = mean(median_primary_err_lo_left_group);
end
end 
%sem 
%sem for individual subject vectors 
primary_err_lo_left_sem = std(primary_err_lo_left_subj)/sqrt(length(primary_err_lo_left_subj));
%primary_err_lo_left_sem = std(median_primary_err_lo_left_group)/sqrt(length(median_primary_err_lo_left_group));


%error bars 
%subj = {'gh','JF','MP','pk','EK','mr','cc'};
primary_err_lo_right_subj = [];
median_primary_err_lo_right_group = [];

for cc =1:length(cond);
for ss = 1:length(subj);
    %filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_lo.mat',subj{ss});
    filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_lo.mat',subj{ss},cond{cc});
    resultsfile = load(filename) %want just ii_results_xx
    primary_err_lo_right_subj = [primary_err_lo_right_subj; resultsfile.ii_results_lo.no_break_right_primary_err_z];
    median_primary_err_lo_right_group = [median_primary_err_lo_right_group; resultsfile.ii_results_lo.median_no_break_right_primary_err_z];
    primary_err_lo_right_group_mean= mean(median_primary_err_lo_right_group);
end
end
%subj = {'gh','JF','MP','pk','EK','mr','cc'};

primary_err_lo_right_sem = std(primary_err_lo_right_subj)/sqrt(length(primary_err_lo_right_subj));


primary_err_hi_left_subj = [];
median_primary_err_hi_left_group= [];

for cc =1:length(cond);
for ss = 1:length(subj);
%filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_hi.mat',subj{ss});
filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
resultsfile = load(filename)
primary_err_hi_left_subj = [primary_err_hi_left_subj; resultsfile.ii_results_hi.no_break_left_primary_err_z_new];
median_primary_err_hi_left_group = [median_primary_err_hi_left_group; resultsfile.ii_results_hi.median_no_break_left_primary_err_z_new];
primary_err_hi_left_group_mean = mean(median_primary_err_hi_left_group);
end 
end
primary_err_hi_left_sem = std(primary_err_hi_left_subj)/sqrt(length(primary_err_hi_left_subj));
%subj = {'gh','JF','MP','pk','EK','mr','cc'};

primary_err_hi_right_subj = [];
median_primary_err_hi_right_group = [];

for cc =1:length(cond);
for ss = 1:length(subj);
filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/%s/ii_results_hi.mat',subj{ss},cond{cc});
%filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/%s/l_ips2/ii_results_hi.mat',subj{ss});
resultsfile = load(filename)
primary_err_hi_right_subj = [primary_err_hi_right_subj; resultsfile.ii_results_hi.no_break_right_primary_err_z_new];
median_primary_err_hi_right_group = [median_primary_err_hi_right_group; resultsfile.ii_results_hi.median_no_break_right_primary_err_z_new];
primary_err_hi_right_group_mean = mean(median_primary_err_hi_right_group);
end 
end 
primary_err_hi_right_sem = std(primary_err_hi_right_subj)/sqrt(length(primary_err_hi_right_subj));

%sem_vect = [primary_err_hi_left_sem primary_err_hi_right_sem primary_err_lo_left_sem primary_err_lo_right_sem];
%barvect = [primary_err_hi_left_group_mean primary_err_hi_right_group_mean primary_err_lo_left_group_mean primary_err_lo_right_group_mean];

%barvect_collapse = [(primary_err_hi_left_group_mean + primary_err_hi_right_group_mean)/2 (primary_err_lo_left_group_mean + primary_err_lo_right_group_mean)/2];



barvect_collapse = [(median(primary_err_hi_left_subj) + median(primary_err_hi_right_subj))/2 (median(primary_err_lo_left_subj) + median(primary_err_lo_right_subj))/2];
barvect_subj = [median(primary_err_hi_left_subj) median(primary_err_lo_left_subj) median(primary_err_hi_right_subj) median(primary_err_lo_right_subj)];
sem_vect = [primary_err_hi_left_sem primary_err_lo_left_sem primary_err_hi_right_sem primary_err_lo_right_sem];



combined_hi = [median_primary_err_hi_left_group; median_primary_err_hi_right_group];



hi_mean = mean(combined_hi)
combined_lo = [median_primary_err_lo_left_group; median_primary_err_lo_right_group];
mean_lo = mean(combined_lo)
both_mean = [hi_mean mean_lo];
%save('C:/Volumes/hyper/experiments/Grace', fileId)
%%
%labels = {'Hi, Left'; 'Hi, Right'; 'Lo, Left'; 'Lo, Right';};
figure(2); 

%plot(1,barvect_subj(1),'o-','color', c1,'markersize',10,'MarkerEdgeColor', c1);
errorbar([barvect_subj(1) barvect_subj(2)], [sem_vect(1) sem_vect(2)],'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1);
%plot(2,barvect_subj(2),'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1);
%errorbar(barvect_subj(2) , sem_vect(2) ,'o--','color', c1,'markersize', 10,'MarkerEdgeColor', c1,'MarkerFaceColor', c1)
hold on;
%plot(3,barvect_subj(3),'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13);
errorbar([barvect_subj(3) barvect_subj(4)] ,[sem_vect(3) sem_vect(4)],'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13);
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel',  {'High Priority', 'Low Priority'})
set(gca, 'fontsize', 14)
ylim([0.75 2.75])
%plot(4,barvect_subj(4),'o--','color', c13,'markersize', 10,'MarkerEdgeColor', c13,'MarkerFaceColor', c13);
%errorbar(barvect_subj(4) , sem_vect(4) ,'o-','color', c13,'markersize',10,'MarkerEdgeColor', c13)

%errorbar(barvect_subj,sem_vect,'color', c1,'markersize',10,'MarkerEdgeColor', c1,'.')
title([sprintf(' %s,',subj{ss}),sprintf(' %s:',cond{cc}) 'Median WM error, primary   '])
ylabel ('Error (DVA)')
set(gca,'FontSize',14);
xlim([0 3])
legend({'Ipsi', 'Contra'}); 

%%
% [vect1 h] = ranksum(median_primary_err_hi_left_group, median_primary_err_lo_left_group)
% [vect2 h] = ranksum(median_primary_err_hi_left_group, median_primary_err_hi_right_group)
% [vect3 h] = ranksum(median_primary_err_hi_right_group, median_primary_err_lo_right_group)
% [vect4 h] = ranksum(median_primary_err_lo_left_group, median_primary_err_lo_right_group)
[vect1 h] = ranksum(primary_err_hi_left_subj, primary_err_lo_left_subj)
[vect2 h] = ranksum(primary_err_hi_left_subj, primary_err_hi_right_subj)
[vect3 h] = ranksum(primary_err_hi_right_subj, primary_err_lo_right_subj)
[vect4 h] = ranksum(primary_err_lo_left_subj, primary_err_lo_right_subj)

sigvect = [vect1 vect3];
sigstar({[1,3],[2,4]},sigvect)




hi_sem = std(combined_hi)/sqrt(length(combined_hi));
lo_sem= std(combined_lo)/sqrt(length(combined_lo)); 
both_sem =[hi_sem lo_sem]; 


subplot(1,2,2)
plot(1,barvect_collapse(1),'o-','MarkerFacecolor',parspec(1,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1]);
hold on;
plot(2,barvect_collapse(2),'o-','MarkerFacecolor',parspec(13,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1]);
%title('Group average WM error, primary saccade dva')
%xlabel ('Priority Hi/Lo')
ylabel ('Error,dva')
set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
legend({'High', 'Low'}); 
ylim([0 2.5])
xlim([0 3])


errorbar(both_mean,both_sem,'.','linewidth',2,'color',c1)
[vect_col h] = ranksum(combined_hi, combined_lo)
sigstar({[1,2]},vect_col)
set(gca, 'Fontsize', 14)

% % sem (error bars mean nothing until ~ 7 subjects)
% for ii = 1:length(subj);
% meansforsem_lo_left = [mean(final_err_lo_left_subj_final{1,ii});]
% end 


% 
% meansforsem_lo_left = [mean(final_err_lo_left_subj_final{1,1}) mean(final_err_lo_left_subj_final{1,2}) mean(final_err_lo_left_subj_final{1,3})] ;
% sem_lo_left = mean(meansforsem_lo_left)/sqrt(size(meansforsem_lo_left,2));
% 
% meansforsem_lo_right = [mean(thingtoplot_lo_right_subj_final{1,1}) mean(thingtoplot_lo_right_subj_final{1,2}) mean(thingtoplot_lo_right_subj_final{1,3})] ;
% sem_lo_right = mean(meansforsem_lo_right)/sqrt(size(meansforsem_lo_right,2));
% 
% meansforsem_hi_left = [mean(thingtoplot_hi_left_subj_final{1,1}) mean(thingtoplot_hi_left_subj_final{1,2}) mean(thingtoplot_hi_left_subj_final{1,3})] ;
% sem_hi_left = mean(meansforsem_hi_left)/sqrt(size(meansforsem_hi_left,2));
% 
% meansforsem_hi_right = [mean(thingtoplot_hi_right_subj_final{1,1}) mean(thingtoplot_hi_right_subj_final{1,2}) mean(thingtoplot_hi_right_subj_final{1,3})] ;
% sem_hi_right = mean(meansforsem_hi_right)/sqrt(size(meansforsem_hi_right,2));

sem_vect = [sem_lo_left sem_lo_right sem_hi_left sem_hi_right];
%% FINAL EYE POSITION
%need to be able to access means of each subj 
% subj = {'gh','JF','MP','pk','EK','mr','cc'};
subj = {'gh','MP','pk','EK','mr','cc'};
final_err_lo_left_subj = {};
median_final_err_lo_left_group = [];

for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_lo.mat',subj{ss});
    resultsfile = load(filename)
    final_err_lo_left_subj = [resultsfile.ii_results_lo.no_break_left_final_err_z];
    median_final_err_lo_left_group = [median_final_err_lo_left_group; resultsfile.ii_results_lo.median_no_break_left_final_err_z_new];
    final_err_lo_left_group_mean = mean(median_final_err_lo_left_group);
end

final_err_lo_left_sem = std(median_final_err_lo_left_group)/sqrt(length(median_final_err_lo_left_group));



%sem 

% meansforsem = [mean(final_err_lo_left_subj_final{1,1}) mean(final_err_lo_left_subj_final{1,2})];
% sem = mean(meansforsem)/sqrt(size(meansforsem,2)); 

%error bars 

subj = {'gh','MP','pk','EK','mr','cc'};
final_err_lo_right_subj = {};
median_final_err_lo_right_group = [];

for ss = 1:length(subj);
    filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_lo.mat',subj{ss});
    resultsfile = load(filename) %want just ii_results_xx
    final_err_lo_right_subj = [resultsfile.ii_results_lo.no_break_right_final_err_z];
    median_final_err_lo_right_group = [median_final_err_lo_right_group; resultsfile.ii_results_lo.median_no_break_right_final_err_z_new];
    final_err_lo_right_group_mean = mean(median_final_err_lo_right_group);
end
final_err_lo_right_sem = std(median_final_err_lo_right_group)/sqrt(length(median_final_err_lo_right_group));
subj = {'gh','MP','pk','EK','mr','cc'};
final_err_hi_left_subj = {};
median_final_err_hi_left_group = [];

for ss = 1:length(subj);
filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_hi.mat',subj{ss});
resultsfile = load(filename)
final_err_hi_left_subj = [resultsfile.ii_results_hi.median_no_break_left_final_err_z];
median_final_err_hi_left_group = [median_final_err_hi_left_group; resultsfile.ii_results_hi.median_no_break_left_final_err_z_new];
final_err_hi_left_group_mean = mean(median_final_err_hi_left_group);
end 


final_err_hi_left_sem = std(median_final_err_hi_left_group)/sqrt(length(median_final_err_hi_left_group));
subj = {'gh','MP','pk','EK','mr','cc'};
final_err_hi_right_subj = {};
median_final_err_hi_right_group_final = [];

for ss = 1:length(subj);
filename = sprintf('/Volumes/hyper/experiments/Grace/DATA/%s/ii_results_hi.mat',subj{ss});
resultsfile = load(filename)
final_err_hi_right_subj = [resultsfile.ii_results_hi.median_no_break_right_final_err_z];
median_final_err_hi_right_group = [median_final_err_hi_right_group_final; resultsfile.ii_results_hi.median_no_break_right_final_err_z_new];
final_err_hi_right_group_mean = mean(median_final_err_hi_right_group);
end 

combined_hi_fin = [median_final_err_hi_left_group; median_final_err_hi_right_group];
hi_mean_fin = mean(combined_hi_fin)
combined_lo_fin = [median_final_err_lo_left_group; median_final_err_lo_right_group];
mean_lo_fin = mean(combined_lo_fin)
both_mean_fin = [hi_mean_fin mean_lo_fin];

hi_sem_fin = std(combined_hi_fin)/sqrt(length(combined_hi_fin));
lo_sem_fin= std(combined_lo_fin)/sqrt(length(combined_lo_fin)); 

both_sem_fin = [hi_sem_fin lo_sem_fin];
errorbar(both_mean_fin,both_sem_fin,'.','linewidth',2)
[vect_col h] = ranksum(combined_hi_fin, combined_lo_fin)
sigstar({[1,2]},vect_col)
set(gca, 'Fontsize', 14)

%final_err_hi_right_sem = std(median_final_err_hi_right_group)/sqrt(length(median_final_err_hi_right_group));
barvect_final = [final_err_hi_left_group_mean final_err_hi_right_group_mean final_err_lo_left_group_mean final_err_lo_right_group_mean];
barvect_collapse_final = [(final_err_hi_left_group_mean + final_err_hi_right_group_mean)/2 (final_err_lo_left_group_mean + final_err_lo_right_group_mean)/2];

%%


%labels = {'Hi, Left'; 'Hi, Right'; 'Lo, Left'; 'Lo, Right';};

figure(2); 
subplot(1,2,1)
bar(1,barvect_final(1),'Facecolor',parspec(1,:,:),'EdgeColor',[1 1 1]);
hold on;
bar(2,barvect_final(2),'Facecolor',parspec(3,:,:),'EdgeColor',[1 1 1]);
bar(3,barvect_final(3),'Facecolor',parspec(13,:,:),'EdgeColor',[1 1 1]);
bar(4,barvect_final(4),'Facecolor',parspec(16,:,:),'EdgeColor',[1 1 1]);
hold on;
errorbar(barvect,sem_vect,'.')
title('Group average WM error, final saccade dva')
%xlabel ('Priority Hi/Lo, L/R')
ylabel ('Error,dva')
set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
legend({'High, Left', 'High, Right', 'Low, Left', 'Low, right',})
%labels = {'Hi, Left'; 'Hi, Right'; 'Lo, Left'; 'Lo, Right';};
subplot(1,2,2)
plot(1, barvect_collapse_final(1),'o-','MarkerFacecolor',parspec(1,:,:),'MarkerSize',15,'MarkerEdgeColor',[1 1 1])
hold on;
plot(2, barvect_collapse_final(2),'o-','MarkerFacecolor',parspec(13,:,:),'MarkerSize',15,'MarkerEdgeColor',[1 1 1])
errorbar(barvect_collapse_final,both_sem_fin,'.')
title('Group average WM error, final saccade dva')
%xlabel ('Priority Hi/Lo')
ylabel ('Error,dva')
set(gca,'FontSize',14);
set(gca,'XTickLabel', []);
legend({'High', 'Low'})
[vect_col h] = ranksum(combined_hi_fin, combined_lo_fin)
sigstar({[1,2]},vect_col)
set(gca, 'Fontsize', 14)
ylim([0 2.5])
xlim([0 3])



% sem (error bars mean nothing until ~ 7 subjects)

meansforsem_lo_left = [mean(final_err_lo_left_subj_final{1,1}) mean(final_err_lo_left_subj_final{1,2}) mean(final_err_lo_left_subj_final{1,3})] ;
sem_lo_left = mean(meansforsem_lo_left)/sqrt(size(meansforsem_lo_left,2));

meansforsem_lo_right = [mean(thingtoplot_lo_right_subj_final{1,1}) mean(thingtoplot_lo_right_subj_final{1,2}) mean(thingtoplot_lo_right_subj_final{1,3})] ;
sem_lo_right = mean(meansforsem_lo_right)/sqrt(size(meansforsem_lo_right,2));

meansforsem_hi_left = [mean(thingtoplot_hi_left_subj_final{1,1}) mean(thingtoplot_hi_left_subj_final{1,2}) mean(thingtoplot_hi_left_subj_final{1,3})] ;
sem_hi_left = mean(meansforsem_hi_left)/sqrt(size(meansforsem_hi_left,2));

meansforsem_hi_right = [mean(thingtoplot_hi_right_subj_final{1,1}) mean(thingtoplot_hi_right_subj_final{1,2}) mean(thingtoplot_hi_right_subj_final{1,3})] ;
sem_hi_right = mean(meansforsem_hi_right)/sqrt(size(meansforsem_hi_right,2));

sem_vect = [sem_lo_left sem_lo_right sem_hi_left sem_hi_right];



%% first year paper fig 
%get both primary and final collapse into one fig 

subplot(1,2,1)
plot(1,barvect_collapse(1),'o-','MarkerFacecolor',parspec(1,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1]);
hold on;
plot(2,barvect_collapse(2),'o-','MarkerFacecolor',parspec(13,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1]);
%errorbar(barvect,sem_vect,'.')
%title('Group average WM error, primary saccade dva')
%xlabel ('Priority Hi/Lo')
ylabel ('Error,dva')
set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
legend({'High', 'Low'}); 
ylim([0 2.5])
%xlim([0 5])
hi_sem = std(combined_hi)/sqrt(length(combined_hi));
lo_sem= std(combined_lo)/sqrt(length(combined_lo)); 

both_sem =[hi_sem lo_sem]; 

errorbar(both_mean,both_sem,'.','linewidth',2)
[vect_col h] = ranksum(combined_hi, combined_lo)
sigstar({[1,2]},vect_col)
set(gca, 'Fontsize', 14)


subplot(1,2,2)
plot(1, barvect_collapse_final(1),'o-','MarkerFacecolor',parspec(1,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1])
hold on;
plot(2, barvect_collapse_final(2),'o-','MarkerFacecolor',parspec(13,:,:),'MarkerSize',12,'MarkerEdgeColor',[1 1 1])
errorbar(barvect_collapse_final,both_sem_fin,'.','linewidth',2) 
title('Group average WM error, final saccade dva')
ylabel ('Error,dva')
set(gca,'FontSize',14);
set(gca,'XTickLabel', []);
legend({'High', 'Low'})
[vect_col h] = ranksum(combined_hi_fin, combined_lo_fin)
sigstar({[1,2]},vect_col)
set(gca, 'Fontsize', 14)
ylim([0 2.5])
xlim([0 3])


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



    
    



