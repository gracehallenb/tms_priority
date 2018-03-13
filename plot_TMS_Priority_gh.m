% based on plot_topoTMS_tcs.m by TCS, 2018 

cd ('/Volumes/hyper/experiments/Grace/TMS_Priority')
data_dir = '/Volumes/hyper/experiments/Grace/TMS_Priority';

subj = {'subj01','subj02', 'subj03','subj04'}; 
cond = {'noTMS','l_spcs','l_ips2'};

contra = [NaN 2 2]; % 1 = LEFT, 2 = RIGHT, NaN = no sorting

nblank = length(subj)*length(cond)*396;

all_subj = nan(nblank,1);
all_cond = nan(nblank,1);
all_hemi = nan(nblank,1); 
all_hemipos =nan(nblank,1); % probably redundant...
all_run  = nan(nblank,1);
all_pri =  nan(nblank,1);
all_targ = nan(nblank,2); %target coords

all_sacc_p = nan(nblank,2); % primary, final saccade points (from ii_stats)
all_sacc_f = nan(nblank,2);

all_err_p = nan(nblank,1); % primary, final error (from ii_stats)
all_err_f = nan(nblank,1);
all_srt   = nan(nblank,1);

all_compliance = nan(nblank,3); % will put ii_stats.trial_compliance here (no idea what's in here, but 0's are bad)

startidx = 1;


%% load everything
for ss = 1:length(subj)
    for cc = 1:length(cond)
        
        % first, load ii_stats
        statsfn = sprintf('%s/%s/%s/ii_stats.mat',data_dir,subj{ss},cond{cc});
        fprintf('Loading ii_stats: %s\n',statsfn);
        load(statsfn);
    
        
        % use # of runs in ii_stats to determine # of runs to load from
        % DATA
        tmp_startidx = startidx;
        for rr = 1:length(ii_stats)
            
            
   
        myf = dir(sprintf('%s/%s/%s/TASK/',data_dir,subj{ss},cond{cc}));
        myf = myf(~[myf.isdir]);
        this_run = load(sprintf('%s/%s/%s/TASK/%s',data_dir,subj{ss},cond{cc},myf(rr).name)); %myscreen, task, tasktiming, stimulus info
        
        % load in target positions
        thisidx = tmp_startidx:(tmp_startidx+size(ii_stats(rr).srt_cursel,1)-1);
        all_targ(thisidx,:) = [this_run.task.targCoords{1}(:,1) this_run.task.targCoords{1}(:,2)]; % want X and Y here targCoords from QUERIED TARG here
        all_run(thisidx) = rr;
        all_pri(thisidx,:) = [this_run.task.conditionAndQueriedTarget];
        all_hemipos(thisidx,:) = [this_run.task.hemiPos];
        tmp_startidx = thisidx(end)+1;
        
        clear this_run thisidx;
        end
        clear tmp_startidx;
        
        thisidx = startidx:(startidx+length(ii_stats)*length(ii_stats(1).srt)-1);    
        
        all_subj(thisidx) = ss;
        all_cond(thisidx) = cc;
        
        all_srt(thisidx)        = cat_field(ii_stats,'srt');
        all_err_p(thisidx)      = cat_field(ii_stats,'primary_err_z');
        all_err_f(thisidx)      = cat_field(ii_stats,'final_err_z');
        
        all_sacc_p(thisidx,:)   = [cat_field(ii_stats,'primary_x') cat_field(ii_stats,'primary_y')];
        all_sacc_f(thisidx,:)   = [cat_field(ii_stats,'final_x')   cat_field(ii_stats,'final_y')];
        
        all_compliance(thisidx,:) = cat_field(ii_stats,'trial_compliance');
        
        
        
        startidx = thisidx(end)+1;
        
        % and clean up
        clear ii_stats statsfn thisidx;
        
    end
end

good_idx = ~isnan(all_subj);
all_subj = all_subj(good_idx);
all_cond = all_cond(good_idx);
all_run  = all_run(good_idx);
all_srt = all_srt(good_idx);
all_err_p = all_err_p(good_idx);
all_err_f = all_err_f(good_idx);
all_sacc_p = all_sacc_p(good_idx,:);
all_sacc_f = all_sacc_f(good_idx,:);
all_compliance = all_compliance(good_idx,:);
all_pri = all_pri(good_idx); 
all_hemipos =all_hemipos(good_idx); 

all_targ = all_targ(good_idx,:); 
all_hemi = 1 + (all_targ(:,1)>0); % 1 = left, 2 = right

% NOTE: need to exclude trials with srt > 900, and sum(trial_compliance,2) < 3
excl_trial = all_srt > 900 | all_srt < 100 | all_err_p > 5 | (sum(all_compliance,2)~=3); % 1 for trials to drop
%excl_trial = all_srt > 900 | all_compliance(:,1)~=1; % 1 for trials to drop

cu = unique(all_cond); %1, 2,3 ...
cond_colors = lines(length(cu));

%% plot mean error per hemi per TMS condition

% first, let's just plot L/R err for primary, final saccade for each
% condition

subj_lege = {'subj01','subj02','subj03','subj04'}; 
subj_col = parula(length(subj));
figure; err_ax = []; srt_ax = [];
priu = [31:32];
hu = unique(all_hemi);

priority_id ={'hi','lo'};
hemifield = {'L', 'R'};

for ii =1:length(hu);
for cc = 1:length(cu)
    % subj x hemi
    
    tmp_err_p = nan(length(subj),2);
    tmp_err_f = nan(length(subj),2);
    tmp_srt   = nan(length(subj),2);
      
    for ss = 1:length(subj);
        
        %break it down based on priority
            for pp = 1:length(priu); %priority, 31 == high, 32 == low
                
                thisidx = all_subj==ss & all_cond==cc & all_hemi==ii & all_pri== priu(pp)  & ~excl_trial;
                
                tmp_err_p(ss,pp) = mean(all_err_p(thisidx));
                tmp_err_p_sem(ss,pp) = std(all_err_p(thisidx))/sqrt(length(all_err_p(thisidx)));
                tmp_err_p_var(ss,pp) = var(all_err_p(thisidx));
                
                tmp_err_f(ss,pp) = mean(all_err_f(thisidx));
                tmp_err_f_sem(ss,pp) = std(all_err_f(thisidx))/sqrt(length(all_err_f(thisidx)));
                tmp_err_f_var(ss,pp) = var(all_err_f(thisidx));
                tmp_srt(ss,pp) = mean(all_srt(thisidx));

            end 
            
        
    end
    

    
    for gg = 1:length(subj);
        figure(1)
        hold on;
        subplot(3,3,cc)
        this_x = [1 2] + (ii-1)*2;
        plot(this_x, tmp_err_p(gg,:).', '-','Color',subj_col(gg,:));
        if gg  == length(subj)
            plot(this_x, mean(tmp_err_p),'o-','MarkerSize',8,'Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'LineWidth',2);
        else
        end
        
        if gg  == length(subj) && cc ==3
            legend(subj_lege, 'mean')
        else
        end
        title(cond{cc});
        if cc == 1
            ylabel('Primary error (\circ)');
        end
        set(gca, 'fontsize',14)
        set(gca,'XTick',[1 2 3 4]);
        set(gca,'XTickLabel',{'ipsi,hi','ipsi,lo','contra,hi','contra, lo'});
      
        
    end
  
    
    
    for vv  = 1:length(subj)
        subplot(3,3,cc+length(cu))
        this_x = [1 2] + (ii-1)*2;
        plot(this_x, tmp_err_f(vv,:).', '-','Color',subj_col(vv,:));
        if vv  == length(subj)
            plot(this_x, mean(tmp_err_f),'o-','MarkerSize',8,'Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'LineWidth',2);
        else
        end
        
        set(gca, 'fontsize',14)
        set(gca,'XTick',[1 2 3 4]);
        set(gca,'XTickLabel',{'ipsi,hi','ipsi,lo','contra,hi','contra, lo'});
        if vv  == length(subj) & cc ==3
            legend(subj_lege, 'mean')
        else
        end 
        title(cond{cc});
        if cc == 1
            ylabel('Final error (\circ)');
        end
        
        hold on;
        ylim([0.5 2.5])
        xlim([0 5])
    end
    
    for kk =1:length(subj)
        subplot(3,3,cc+length(cu)*2)
        plot(this_x, tmp_srt(kk,:).', '-','Color',subj_col(kk,:));
        hold on;
        
        if kk == length(subj)
        plot(this_x, mean(tmp_srt),'o-','MarkerSize',8,'Color',cond_colors(cc,:),'MarkerFaceColor',cond_colors(cc,:),'LineWidth',2);
        else 
        end 
        set(gca,'XTick',[0 1 2 3 4 5]);
        set(gca,'XTickLabel',{'ipsi,hi','ipsi,lo','contra,hi','contra, lo'});
        if kk  == length(subj) & cc ==3
            legend(subj_lege, 'mean')
        else
        end 
         if cc == 1
            ylabel('SRT');
        end
        set(gca, 'fontsize',14)
        set(gca,'XTick',[1 2 3 4]);
        xlim([0 5])
    end
    
end 


end



%% scatterplots like Fig. 2 (from Mackey & Curtis 2017)
figure;
cond_scatter_y = [2 3]; % different y axes
cond_scatter_x = [1 1]; % and consistent x axis
for cc = 1:length(cond_scatter_y)
    
    x_err_p = nan(length(subj),2);
    x_err_f = nan(length(subj),2);
    x_srt   = nan(length(subj),2);

    y_err_p = nan(length(subj),2);
    y_err_f = nan(length(subj),2);
    y_srt   = nan(length(subj),2);

    
    for ss = 1:length(subj)
        for ii = 1:2 % L, R
          
            thisidx = all_subj==ss & all_cond==cond_scatter_y(cc) & all_hemi==ii & ~excl_trial;
            y_err_p(ss,ii) = mean(all_err_p(thisidx));
            y_err_f(ss,ii) = mean(all_err_f(thisidx));
            y_srt(ss,ii) = mean(all_srt(thisidx));
            
            thisidx = all_subj==ss & all_cond==cond_scatter_x(cc) & all_hemi==ii & ~excl_trial;
            x_err_p(ss,ii) = mean(all_err_p(thisidx));
            x_err_f(ss,ii) = mean(all_err_f(thisidx));
            x_srt(ss,ii)   = mean(all_srt(thisidx));

        end
    end

    
    
    % primary
    subplot(length(cond_scatter_y),2,1+(cc-1)*length(cond_scatter_y)); hold on;
    plot(x_err_p(:,1),y_err_p(:,1),'o','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',[1 1 1],'Color',cond_colors(cond_scatter_y(cc),:));
    plot(x_err_p(:,2),y_err_p(:,2),'o','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',cond_colors(cond_scatter_y(cc),:),'Color',cond_colors(cond_scatter_y(cc),:));
    plot([0 2.5],[0 2.5],'r--');
    
    axis  equal square;
    xlim([0.8 2.2]);ylim([0.8 2.2]);
    
    xlabel(sprintf('Err: %s',cond{cond_scatter_x(cc)}));
    ylabel(sprintf('Err: %s',cond{cond_scatter_y(cc)}));
    title('Primary error (\circ)');
    set(gca,'XTick',[1 1.5 2],'YTick',[1 1.5 2]);
    
    % final
    subplot(length(cond_scatter_y),2,2+(cc-1)*length(cond_scatter_y)); hold on;
    plot(x_err_f(:,1),y_err_f(:,1),'o','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',[1 1 1],'Color',cond_colors(cond_scatter_y(cc),:));
    plot(x_err_f(:,2),y_err_f(:,2),'o','LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',cond_colors(cond_scatter_y(cc),:),'Color',cond_colors(cond_scatter_y(cc),:));
    plot([0 2.5],[0 2.5],'r--');
    axis  equal square;
    xlim([0.6 1.4]);ylim([0.6 1.4]);
    
    xlabel(sprintf('Err: %s',cond{cond_scatter_x(cc)}));
    ylabel(sprintf('Err: %s',cond{cond_scatter_y(cc)}));
    title('Final error (\circ)');
    set(gca,'XTick',[0.6 1 1.4],'YTick',[0.6 1 1.4]);

end
