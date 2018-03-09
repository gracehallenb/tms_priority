% Adaptation of MGS STATS(2015). Accounts for different saccadic conditions
% by incorporating the function cat_structs (TCS,2017)
%%%%%%%%%%%%%%%%
% WHAT WE WANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% SRT (DONE)
% ERROR (DONE)
% GAIN (DONE)
% VALID TRIALS MATRIX (DONE)
% PERCENT TRIALS BROKE FIXATION (0 = broke, 1 = didnt break)
% BREAK UP STATS IN CATEGORIES: ALL, BROKE FIXATION, DIDN'T BREAK FIXATION
% PLOT DISTRIBUTIONS
% SPLIT DATA BY HEMIFIELD
% SPLIT DATA BY QUADRANT

condition ={'lo'};

for cc =1:length(condition)
%%% which cond? hi or lo? %
cond = sprintf('%s',condition{cc});

%% load ii_stats file 

load('/Volumes/hyper/experiments/Grace/TMS_Priority/subj04/l_ips2/ii_stats.mat');

%runs = {'run21','run23','run24','run25','run26','run02','run03','run04','run05'}; %ek spcs
%runs = {'run01','run02','run03','run04','run05','run07','run08','run09'};
%runs = {'run01','run02','run03','run04','run05','run06','run07','run08','run09','run10'};
%runs = {'run07','run08','run09','run10','run11'};
%runs = {'run02','run03','run04','run05','run06','run07','run08'};
%runs = {'run01','run02','run03','run04','run05','run06','run07','run08','run09'}; %ekpi
%runs ={'run01','run02','run03','run04','run05','run06','run07','run09','run10'}; %ek sham
%runs = {'run01','run02','run03','run04','run05','run06'};mr pilot
%runs = {'run02','run03','run04','run05','run06','run07','run08','run09','run10'};
%runs ={'run01','run02','run03','run04','run05','run06','run07','run09','run10'}; %mr 
%runs = {'run01','run02','run03','run04','run05','run08','run09','run10','run11'}; %ab sham
%runs = {'run01','run03','run04','run05','run12','run13','run14','run15','run16'}; %ab spcs
%runs = {'run01','run02','run07','run09','run10','run11','run12'} %ab lips2
%runs = {'run01','run02','run03','run04','run05','run06'}; %ab pilot
%runs = {'run11', 'run12', 'run13', 'run14','run15', 'run16', 'run17','run18', 'run19', 'run20'}; %eklips2
%runs = {'run12','run13','run14','run15','run16','run01','run03','run04','run05'};
%runs  = {'run07','run09','run10','run11','run12'};
%runs = {'run07','run09','run10','run11','run12','run01','run02'}
%runs ={'run06','run07','run08','run09','run10'};
%runs ={'run01','run03','run04','run05','run06','run07'};
%runs ={'run01','run02','run03','run04','run05','run06','run07','run08'};
%runs
%={'run01','run02','run03','run04','run05','run06','run07','run08'};
%runs ={'run01','run02','run03','run04','run05','run06','run07','run08','run09','run10'};
%runs = {'run11','run12'};
%runs = {'run01','run02','run03','run04','run05','run06','run07','run08'};
%runs = {'run01','run02','run03','run04','run05'};
%runs = {'run11','run12','run22','run23','run24','run25','run26','run27','run28'}; %mr lips2 combined
%runs = {'run13','run14','run15','run16'};
%runs = {'run06','run07','run08','run09','run10'};
runs= {'run12','run13','run14','run15'};
%runs = {'run26','run27','run28'};%newest mr lips
%runs = {'run22','run23','run24','run25'} % second mr lips
%runs = {'run06','run07','run08','run09','run10'};

% all_runs = run_record;
% runs = all_runs.subj04_ips2_II; 
fileweneed = [];
%task_priority_identity = [];  %go into .mat file which was generated during experiment, take the priority cond on each trial
for jj = 1:length(runs);
filename = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj04/l_ips2/TASK/%s.mat',runs{jj});
fileID = load(filename);
% task_priority_identity = [task_priority_identity; fileID.task.conditionAndQueriedTarget(:,1)];
% end 
if  jj == 1;
    newrow = 1;
    endrow = 36;
else
newrow = ((jj-1).*36)+1
endrow = newrow + 35
end; 
 fileweneed(newrow:endrow,1) = fileID.task.conditionAndQueriedTarget(:,1);
end

% concatenate all condition labels

% concatenate ii_stats array into one variable (all_ii)
all_ii = [];

skip_fields = {'trial_notes', 'srt_sel','raw_x','raw_y','acc_sel','acc_cursel'};

for ii = 1:length(ii_stats);
    
    all_ii = cat_struct(all_ii,ii_stats(ii),skip_fields);
    
end

all_ii = rmfield(all_ii, skip_fields);

%%  separate ii_stats data into hi and lo conditions

if cond == 'hi';
    run_num_hi_ind = find(fileweneed == 31); % hi priority identity, should occur 24x/ run
    numtrials_givenpriority = run_num_hi_ind;
else
    run_num_lo_ind = find(fileweneed == 32); %lo priority identity, should occur 12x /run
    numtrials_givenpriority = run_num_lo_ind;
end

ii_stats = [];
 myf = fieldnames(all_ii);
 
 for xx = 1:length(myf);
    ii_stats.(myf{xx}) = all_ii.(myf{xx})(numtrials_givenpriority); %how to index into hi/lo from myf loop
 end


%%
    

trialnum = length(numtrials_givenpriority); 
num_runs = length(ii_stats); % How many runs did the subject complete?

%%%%%%%%%%%%%%%%%%%%%%%%
% TOSS OUT NOGO TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Throw out nogo trials (SRT > 900ms);

ii_results.num_trials = num_runs * trialnum;
ii_results.num_nogo = 0;
ii_results.num_trash = 0;

for i = 1:num_runs
    ii_stats(i).nogo_mat = ones(trialnum,1);
    ii_stats(i).trash_mat = ones(trialnum,1);
    ii_stats(i).nogo_inds = find(ii_stats(i).srt>=900);
    %ii_stats(i).trash_inds = find(ii_stats(i).trash==0);
    ii_stats(i).num_nogo = trialnum - length(ii_stats(i).nogo_inds); % Get # of nogo trials
    %ii_stats(i).num_trash = trialnum - length(ii_stats(i).trash_inds); % Get # of nogo trials
    %ii_stats(i).trash_mat(ii_stats(i).trash_inds) = 0;
    ii_stats(i).nogo_mat(ii_stats(i).nogo_inds) = 0;
    %ii_stats(i).nogo_mat = ii_stats(i).nogo_mat.*ii_stats(i).trash_mat;
    ii_results.num_nogo = ii_results.num_nogo + ii_stats(i).num_nogo;
    %ii_results.num_trash = ii_results.num_trash + ii_stats(i).num_trash;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE VALIDITY MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

ii_results.num_bfix = 0;

% Find fixation breaks
for i = 1:num_runs
    ii_stats(i).break_mat = ii_stats(i).trial_compliance(:,1,1); % Find fixation breaks
    ii_stats(i).num_fix_breaks = trialnum - sum(ii_stats(i).break_mat); % Get # of fixation break trials
    ii_results.num_bfix = ii_results.num_bfix + ii_stats(i).num_fix_breaks;
end

% Find left hemifield trials
for i = 1:num_runs
    ii_stats(i).left_mat = zeros(trialnum,1);
    ii_stats(i).left_inds = find(ii_stats(i).corrective_x>=0);
    ii_stats(i).left_mat(ii_stats(i).left_inds) = 1;
    ii_stats(i).num_left_trials = trialnum - length(ii_stats(i).left_inds); 
end

% Find right hemifield trials
for i = 1:num_runs
    ii_stats(i).right_mat = zeros(trialnum,1);
    ii_stats(i).right_inds = find(ii_stats(i).corrective_x<0);
    ii_stats(i).right_mat(ii_stats(i).right_inds) = 1; 
    ii_stats(i).num_right_trials = trialnum - length(ii_stats(i).right_inds);
end
%%
% Create NO fixation break matrix
for i = 1:num_runs
    ii_stats(i).no_break_mat = zeros(trialnum,1);
    ii_stats(i).good_inds = find(ii_stats(i).nogo_mat==1);
    ii_stats(i).no_break_inds = find(ii_stats(i).break_mat==1);
    ii_stats(i).no_break_mat(ii_stats(i).no_break_inds) = 1;
    ii_stats(i).no_break_mat = ii_stats(i).no_break_mat .* ii_stats(i).nogo_mat;
    ii_stats(i).no_break_left_mat = ii_stats(i).no_break_mat .* ii_stats(i).left_mat;
    ii_stats(i).no_break_right_mat = ii_stats(i).no_break_mat .* ii_stats(i).right_mat;
end

% Create ONLY fixation break matrix
for i = 1:num_runs
    ii_stats(i).only_break_mat = zeros(trialnum,1);
    ii_stats(i).good_inds = find(ii_stats(i).nogo_mat==1);
    ii_stats(i).only_break_inds = find(ii_stats(i).break_mat==0);
    ii_stats(i).only_break_mat(ii_stats(i).only_break_inds) = 1;
    ii_stats(i).only_break_mat = ii_stats(i).only_break_mat .* ii_stats(i).nogo_mat;
    ii_stats(i).only_break_left_mat = ii_stats(i).only_break_mat .* ii_stats(i).left_mat;
    ii_stats(i).only_break_right_mat = ii_stats(i).only_break_mat .* ii_stats(i).right_mat;
end

%%%%%%%%%%%%%%%%%%%%
% COMBINE ALL RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW STATS (NOGO ONLY) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% All trials
for i = 1:num_runs
    ii_stats(i).good_inds = find(ii_stats(i).nogo_mat==1); % Get good trials
    
    % Error
    ii_stats(i).good_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).good_inds);
    ii_stats(i).good_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).good_inds);
    ii_stats(i).good_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).good_inds);
    
    ii_stats(i).good_final_err_x = ii_stats(i).final_err_x(ii_stats(i).good_inds);
    ii_stats(i).good_final_err_y = ii_stats(i).final_err_y(ii_stats(i).good_inds);
    ii_stats(i).good_final_err_z = ii_stats(i).final_err_z(ii_stats(i).good_inds);
    
    % Gain
    ii_stats(i).good_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).good_inds);
    ii_stats(i).good_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).good_inds);
    ii_stats(i).good_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).good_inds);
    
    ii_stats(i).good_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).good_inds);
    ii_stats(i).good_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).good_inds);
    ii_stats(i).good_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).good_inds);
    
    % SRT
    ii_stats(i).good_srt = ii_stats(i).srt(ii_stats(i).good_inds);
end

% Left trials
for i = 1:num_runs
    ii_stats(i).left_inds = find(ii_stats(i).left_mat==1); % Get good trials
    
    % Error
    ii_stats(i).left_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).left_inds);
    ii_stats(i).left_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).left_inds);
    ii_stats(i).left_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).left_inds);
    
    ii_stats(i).left_final_err_x = ii_stats(i).final_err_x(ii_stats(i).left_inds);
    ii_stats(i).left_final_err_y = ii_stats(i).final_err_y(ii_stats(i).left_inds);
    ii_stats(i).left_final_err_z = ii_stats(i).final_err_z(ii_stats(i).left_inds);
    
    % Gain
    ii_stats(i).left_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).left_inds);
    ii_stats(i).left_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).left_inds);
    ii_stats(i).left_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).left_inds);
    
    ii_stats(i).left_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).left_inds);
    ii_stats(i).left_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).left_inds);
    ii_stats(i).left_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).left_inds);
    
    % SRT
    ii_stats(i).left_srt = ii_stats(i).srt(ii_stats(i).left_inds);
end

% Right trials
for i = 1:num_runs
    ii_stats(i).right_inds = find(ii_stats(i).right_mat==1); % Get good trials
    
    % Error
    ii_stats(i).right_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).right_inds);
    ii_stats(i).right_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).right_inds);
    ii_stats(i).right_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).right_inds);
    
    ii_stats(i).right_final_err_x = ii_stats(i).final_err_x(ii_stats(i).right_inds);
    ii_stats(i).right_final_err_y = ii_stats(i).final_err_y(ii_stats(i).right_inds);
    ii_stats(i).right_final_err_z = ii_stats(i).final_err_z(ii_stats(i).right_inds);
    
    % Gain
    ii_stats(i).right_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).right_inds);
    ii_stats(i).right_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).right_inds);
    ii_stats(i).right_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).right_inds);
    
    ii_stats(i).right_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).right_inds);
    ii_stats(i).right_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).right_inds);
    ii_stats(i).right_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).right_inds);
    
    % SRT
    ii_stats(i).right_srt = ii_stats(i).srt(ii_stats(i).right_inds);
end

% All trials
ii_results.all_primary_err_x = [];
ii_results.all_primary_err_y = [];
ii_results.all_primary_err_z = [];

ii_results.all_final_err_x = [];
ii_results.all_final_err_y = [];
ii_results.all_final_err_z = [];

ii_results.all_primary_gain_x = [];
ii_results.all_primary_gain_y = [];
ii_results.all_primary_gain_z = [];

ii_results.all_final_gain_x = [];
ii_results.all_final_gain_y = [];
ii_results.all_final_gain_z = [];

ii_results.all_srt = [];

% Left trials
ii_results.all_left_primary_err_x = [];
ii_results.all_left_primary_err_y = [];
ii_results.all_left_primary_err_z = [];

ii_results.all_left_final_err_x = [];
ii_results.all_left_final_err_y = [];
ii_results.all_left_final_err_z = [];

ii_results.all_left_primary_gain_x = [];
ii_results.all_left_primary_gain_y = [];
ii_results.all_left_primary_gain_z = [];

ii_results.all_left_final_gain_x = [];
ii_results.all_left_final_gain_y = [];
ii_results.all_left_final_gain_z = [];

ii_results.all_left_srt = [];

% Right trials
ii_results.all_right_primary_err_x = [];
ii_results.all_right_primary_err_y = [];
ii_results.all_right_primary_err_z = [];

ii_results.all_right_final_err_x = [];
ii_results.all_right_final_err_y = [];
ii_results.all_right_final_err_z = [];

ii_results.all_right_primary_gain_x = [];
ii_results.all_right_primary_gain_y = [];
ii_results.all_right_primary_gain_z = [];

ii_results.all_right_final_gain_x = [];
ii_results.all_right_final_gain_y = [];
ii_results.all_right_final_gain_z = [];

ii_results.all_right_srt = [];

% All trials
for j = 1:num_runs
    % Error
    ii_results.all_primary_err_x = [ii_results.all_primary_err_x; ii_stats(j).good_primary_err_x];
    ii_results.all_primary_err_y = [ii_results.all_primary_err_y; ii_stats(j).good_primary_err_y];
    ii_results.all_primary_err_z = [ii_results.all_primary_err_z; ii_stats(j).good_primary_err_z];

    ii_results.all_final_err_x = [ii_results.all_final_err_x; ii_stats(j).good_final_err_x];
    ii_results.all_final_err_y = [ii_results.all_final_err_y; ii_stats(j).good_final_err_y];
    ii_results.all_final_err_z = [ii_results.all_final_err_z; ii_stats(j).good_final_err_z];
    
    % Gain
    ii_results.all_primary_gain_x = [ii_results.all_primary_gain_x; ii_stats(j).good_primary_gain_x];
    ii_results.all_primary_gain_y = [ii_results.all_primary_gain_y; ii_stats(j).good_primary_gain_y];
    ii_results.all_primary_gain_z = [ii_results.all_primary_gain_z; ii_stats(j).good_primary_gain_z];
    
    ii_results.all_final_gain_x = [ii_results.all_final_gain_x; ii_stats(j).good_final_gain_x];
    ii_results.all_final_gain_y = [ii_results.all_final_gain_y; ii_stats(j).good_final_gain_y];
    ii_results.all_final_gain_z = [ii_results.all_final_gain_z; ii_stats(j).good_final_gain_z];
    
    % SRT
    
    ii_results.all_srt = [ii_results.all_srt; ii_stats(j).good_srt];
end

% Left trials
for j = 1:num_runs
    % Error
    ii_results.all_left_primary_err_x = [ii_results.all_left_primary_err_x; ii_stats(j).left_primary_err_x];
    ii_results.all_left_primary_err_y = [ii_results.all_left_primary_err_y; ii_stats(j).left_primary_err_y];
    ii_results.all_left_primary_err_z = [ii_results.all_left_primary_err_z; ii_stats(j).left_primary_err_z];

    ii_results.all_left_final_err_x = [ii_results.all_left_final_err_x; ii_stats(j).left_final_err_x];
    ii_results.all_left_final_err_y = [ii_results.all_left_final_err_y; ii_stats(j).left_final_err_y];
    ii_results.all_left_final_err_z = [ii_results.all_left_final_err_z; ii_stats(j).left_final_err_z];
    
    % Gain
    ii_results.all_left_primary_gain_x = [ii_results.all_left_primary_gain_x; ii_stats(j).left_primary_gain_x];
    ii_results.all_left_primary_gain_y = [ii_results.all_left_primary_gain_y; ii_stats(j).left_primary_gain_y];
    ii_results.all_left_primary_gain_z = [ii_results.all_left_primary_gain_z; ii_stats(j).left_primary_gain_z];
    
    ii_results.all_left_final_gain_x = [ii_results.all_left_final_gain_x; ii_stats(j).left_final_gain_x];
    ii_results.all_left_final_gain_y = [ii_results.all_left_final_gain_y; ii_stats(j).left_final_gain_y];
    ii_results.all_left_final_gain_z = [ii_results.all_left_final_gain_z; ii_stats(j).left_final_gain_z];
    
    % SRT
    
    ii_results.all_left_srt = [ii_results.all_left_srt; ii_stats(j).left_srt];
end

% Right trials
for j = 1:num_runs
    % Error
    ii_results.all_right_primary_err_x = [ii_results.all_right_primary_err_x; ii_stats(j).right_primary_err_x];
    ii_results.all_right_primary_err_y = [ii_results.all_right_primary_err_y; ii_stats(j).right_primary_err_y];
    ii_results.all_right_primary_err_z = [ii_results.all_right_primary_err_z; ii_stats(j).right_primary_err_z];

    ii_results.all_right_final_err_x = [ii_results.all_right_final_err_x; ii_stats(j).right_final_err_x];
    ii_results.all_right_final_err_y = [ii_results.all_right_final_err_y; ii_stats(j).right_final_err_y];
    ii_results.all_right_final_err_z = [ii_results.all_right_final_err_z; ii_stats(j).right_final_err_z];
    
    % Gain
    ii_results.all_right_primary_gain_x = [ii_results.all_right_primary_gain_x; ii_stats(j).right_primary_gain_x];
    ii_results.all_right_primary_gain_y = [ii_results.all_right_primary_gain_y; ii_stats(j).right_primary_gain_y];
    ii_results.all_right_primary_gain_z = [ii_results.all_right_primary_gain_z; ii_stats(j).right_primary_gain_z];
    
    ii_results.all_right_final_gain_x = [ii_results.all_right_final_gain_x; ii_stats(j).right_final_gain_x];
    ii_results.all_right_final_gain_y = [ii_results.all_right_final_gain_y; ii_stats(j).right_final_gain_y];
    ii_results.all_right_final_gain_z = [ii_results.all_right_final_gain_z; ii_stats(j).right_final_gain_z];
    
    % SRT
    
    ii_results.all_right_srt = [ii_results.all_right_srt; ii_stats(j).right_srt];
end

% All trials
ii_results.mean_all_srt = mean(ii_results.all_srt);
ii_results.mean_all_primary_err_z = mean(ii_results.all_primary_err_z);
ii_results.mean_all_final_err_z = mean(ii_results.all_final_err_z);
ii_results.mean_all_primary_gain_z = mean(ii_results.all_primary_gain_z);
ii_results.mean_all_final_gain_z = mean(ii_results.all_final_gain_z);

ii_results.median_all_srt = median(ii_results.all_srt);
ii_results.median_all_primary_err_z = median(ii_results.all_primary_err_z);
ii_results.median_all_final_err_z = median(ii_results.all_final_err_z);
ii_results.median_all_primary_gain_z = median(ii_results.all_primary_gain_z);
ii_results.median_all_final_gain_z = median(ii_results.all_final_gain_z);

ii_results.std_all_srt = std(ii_results.all_srt);
ii_results.std_all_primary_err_z = std(ii_results.all_primary_err_z);
ii_results.std_all_final_err_z = std(ii_results.all_final_err_z);
ii_results.std_all_primary_gain_z = std(ii_results.all_primary_gain_z);
ii_results.std_all_final_gain_z = std(ii_results.all_final_gain_z);

% Left trials
ii_results.mean_all_left_srt = mean(ii_results.all_left_srt);
ii_results.mean_all_left_primary_err_z = mean(ii_results.all_left_primary_err_z);
ii_results.mean_all_left_final_err_z = mean(ii_results.all_left_final_err_z);
ii_results.mean_all_left_primary_gain_z = mean(ii_results.all_left_primary_gain_z);
ii_results.mean_all_left_final_gain_z = mean(ii_results.all_left_final_gain_z);

ii_results.median_all_left_srt = median(ii_results.all_left_srt);
ii_results.median_all_left_primary_err_z = median(ii_results.all_left_primary_err_z);
ii_results.median_all_left_final_err_z = median(ii_results.all_left_final_err_z);
ii_results.median_all_left_primary_gain_z = median(ii_results.all_left_primary_gain_z);
ii_results.median_all_left_final_gain_z = median(ii_results.all_left_final_gain_z);

ii_results.std_all_left_srt = std(ii_results.all_left_srt);
ii_results.std_all_left_primary_err_z = std(ii_results.all_left_primary_err_z);
ii_results.std_all_left_final_err_z = std(ii_results.all_left_final_err_z);
ii_results.std_all_left_primary_gain_z = std(ii_results.all_left_primary_gain_z);
ii_results.std_all_left_final_gain_z = std(ii_results.all_left_final_gain_z);

% Right trials
ii_results.mean_all_right_srt = mean(ii_results.all_right_srt);
ii_results.mean_all_right_primary_err_z = mean(ii_results.all_right_primary_err_z);
ii_results.mean_all_right_final_err_z = mean(ii_results.all_right_final_err_z);
ii_results.mean_all_right_primary_gain_z = mean(ii_results.all_right_primary_gain_z);
ii_results.mean_all_right_final_gain_z = mean(ii_results.all_right_final_gain_z);

ii_results.median_all_right_srt = median(ii_results.all_right_srt);
ii_results.median_all_right_primary_err_z = median(ii_results.all_right_primary_err_z);
ii_results.median_all_right_final_err_z = median(ii_results.all_right_final_err_z);
ii_results.median_all_right_primary_gain_z = median(ii_results.all_right_primary_gain_z);
ii_results.median_all_right_final_gain_z = median(ii_results.all_right_final_gain_z);

ii_results.std_all_right_srt = std(ii_results.all_right_srt);
ii_results.std_all_right_primary_err_z = std(ii_results.all_right_primary_err_z);
ii_results.std_all_right_final_err_z = std(ii_results.all_right_final_err_z);
ii_results.std_all_right_primary_gain_z = std(ii_results.all_right_primary_gain_z);
ii_results.std_all_right_final_gain_z = std(ii_results.all_right_final_gain_z);


%%%%%%%%%%%%%%%%%%%%%%
% NO FIXATION BREAKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

% All trials
for i = 1:num_runs
    ii_stats(i).no_break_inds = find(ii_stats(i).no_break_mat==1); % Get good trials
    
    % Error
    ii_stats(i).no_break_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_inds);
    
    ii_stats(i).no_break_final_err_x = ii_stats(i).final_err_x(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_final_err_y = ii_stats(i).final_err_y(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_inds);
    
    % Gain
    ii_stats(i).no_break_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_inds);
    
    ii_stats(i).no_break_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).no_break_inds);
    ii_stats(i).no_break_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_inds);
    
    % SRT
    ii_stats(i).no_break_srt = ii_stats(i).srt(ii_stats(i).no_break_inds);
    
%     % MS
%     ii_stats(i).no_break_ms_duration = ii_stats(i).ms_duration(ii_stats(i).no_break_inds);
%     ii_stats(i).no_break_ms_peak_velocity = ii_stats(i).ms_peak_velocity(ii_stats(i).no_break_inds);
   %ii_stats(i).no_break_primary_theta = ii_stats(i).primary_theta(ii_stats(i).no_break_inds);
end

% Left trials
for i = 1:num_runs
    ii_stats(i).no_break_left_inds = find(ii_stats(i).no_break_left_mat==1); % Get good trials
    
    % Error
    ii_stats(i).no_break_left_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_left_inds);
    
    ii_stats(i).no_break_left_final_err_x = ii_stats(i).final_err_x(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_final_err_y = ii_stats(i).final_err_y(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_left_inds);
    
    % Gain
    ii_stats(i).no_break_left_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_left_inds);
    
    ii_stats(i).no_break_left_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).no_break_left_inds);
    ii_stats(i).no_break_left_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_left_inds);
    
    % SRT
    ii_stats(i).no_break_left_srt = ii_stats(i).srt(ii_stats(i).no_break_left_inds);
    
    %theta
     %ii_stats(i).no_break_left_primary_theta = ii_stats(i).primary_theta(ii_stats(i).no_break_left_inds);
end

% Right trials
for i = 1:num_runs
    ii_stats(i).no_break_right_inds = find(ii_stats(i).no_break_right_mat==1); % Get good trials
    
    % Error
    ii_stats(i).no_break_right_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_right_inds);
    
    ii_stats(i).no_break_right_final_err_x = ii_stats(i).final_err_x(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_final_err_y = ii_stats(i).final_err_y(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_right_inds);
    
    % Gain
    ii_stats(i).no_break_right_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_right_inds);
    
    ii_stats(i).no_break_right_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).no_break_right_inds);
    ii_stats(i).no_break_right_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_right_inds);
    
    % SRT
    ii_stats(i).no_break_right_srt = ii_stats(i).srt(ii_stats(i).no_break_right_inds);
     
    %theta
    %ii_stats(i).no_break_right_primary_theta = ii_stats(i).primary_theta(ii_stats(i).no_break_right_inds);
end

% All trials
ii_results.no_break_primary_err_x = [];
ii_results.no_break_primary_err_y = [];
ii_results.no_break_primary_err_z = [];

ii_results.no_break_final_err_x = [];
ii_results.no_break_final_err_y = [];
ii_results.no_break_final_err_z = [];

ii_results.no_break_primary_gain_x = [];
ii_results.no_break_primary_gain_y = [];
ii_results.no_break_primary_gain_z = [];

ii_results.no_break_final_gain_x = [];
ii_results.no_break_final_gain_y = [];
ii_results.no_break_final_gain_z = [];

ii_results.no_break_srt = [];

ii_results.no_break_ms_duration = [];
ii_results.no_break_ms_peak_velocity = [];
%ii_results.no_break_primary_theta = [];

% Left trials
ii_results.no_break_left_primary_err_x = [];
ii_results.no_break_left_primary_err_y = [];
ii_results.no_break_left_primary_err_z = [];

ii_results.no_break_left_final_err_x = [];
ii_results.no_break_left_final_err_y = [];
ii_results.no_break_left_final_err_z = [];

ii_results.no_break_left_primary_gain_x = [];
ii_results.no_break_left_primary_gain_y = [];
ii_results.no_break_left_primary_gain_z = [];

ii_results.no_break_left_final_gain_x = [];
ii_results.no_break_left_final_gain_y = [];
ii_results.no_break_left_final_gain_z = [];

ii_results.no_break_left_srt = [];

% Right trials
ii_results.no_break_right_primary_err_x = [];
ii_results.no_break_right_primary_err_y = [];
ii_results.no_break_right_primary_err_z = [];

ii_results.no_break_right_final_err_x = [];
ii_results.no_break_right_final_err_y = [];
ii_results.no_break_right_final_err_z = [];

ii_results.no_break_right_primary_gain_x = [];
ii_results.no_break_right_primary_gain_y = [];
ii_results.no_break_right_primary_gain_z = [];

ii_results.no_break_right_final_gain_x = [];
ii_results.no_break_right_final_gain_y = [];
ii_results.no_break_right_final_gain_z = [];

ii_results.no_break_right_srt = [];

% All trials
for j = 1:num_runs
    % Error
    ii_results.no_break_primary_err_x = [ii_results.no_break_primary_err_x; ii_stats(j).no_break_primary_err_x];
    ii_results.no_break_primary_err_y = [ii_results.no_break_primary_err_y; ii_stats(j).no_break_primary_err_y];
    ii_results.no_break_primary_err_z = [ii_results.no_break_primary_err_z; ii_stats(j).no_break_primary_err_z];

    ii_results.no_break_final_err_x = [ii_results.no_break_final_err_x; ii_stats(j).no_break_final_err_x];
    ii_results.no_break_final_err_y = [ii_results.no_break_final_err_y; ii_stats(j).no_break_final_err_y];
    ii_results.no_break_final_err_z = [ii_results.no_break_final_err_z; ii_stats(j).no_break_final_err_z];
    
    % Gain
    ii_results.no_break_primary_gain_x = [ii_results.no_break_primary_gain_x; ii_stats(j).no_break_primary_gain_x];
    ii_results.no_break_primary_gain_y = [ii_results.no_break_primary_gain_y; ii_stats(j).no_break_primary_gain_y];
    ii_results.no_break_primary_gain_z = [ii_results.no_break_primary_gain_z; ii_stats(j).no_break_primary_gain_z];
    
    ii_results.no_break_final_gain_x = [ii_results.no_break_final_gain_x; ii_stats(j).no_break_final_gain_x];
    ii_results.no_break_final_gain_y = [ii_results.no_break_final_gain_y; ii_stats(j).no_break_final_gain_y];
    ii_results.no_break_final_gain_z = [ii_results.no_break_final_gain_z; ii_stats(j).no_break_final_gain_z];
    
    % SRT    
    ii_results.no_break_srt = [ii_results.no_break_srt; ii_stats(j).no_break_srt];
    
%     % MS    
%     ii_results.no_break_ms_duration = [ii_results.no_break_ms_duration; ii_stats(j).no_break_ms_duration];
%     ii_results.no_break_ms_peak_velocity = [ii_results.no_break_ms_peak_velocity; ii_stats(j).no_break_ms_peak_velocity];
%theta 
    %ii_results.no_break_primary_theta = [ii_results.no_break_primary_theta; ii_stats(j).no_break_primary_theta];
end

% Left trials
for j = 1:num_runs
    % Error
    ii_results.no_break_left_primary_err_x = [ii_results.no_break_left_primary_err_x; ii_stats(j).no_break_left_primary_err_x];
    ii_results.no_break_left_primary_err_y = [ii_results.no_break_left_primary_err_y; ii_stats(j).no_break_left_primary_err_y];
    ii_results.no_break_left_primary_err_z = [ii_results.no_break_left_primary_err_z; ii_stats(j).no_break_left_primary_err_z];

    ii_results.no_break_left_final_err_x = [ii_results.no_break_left_final_err_x; ii_stats(j).no_break_left_final_err_x];
    ii_results.no_break_left_final_err_y = [ii_results.no_break_left_final_err_y; ii_stats(j).no_break_left_final_err_y];
    ii_results.no_break_left_final_err_z = [ii_results.no_break_left_final_err_z; ii_stats(j).no_break_left_final_err_z];
    
    % Gain
    ii_results.no_break_left_primary_gain_x = [ii_results.no_break_left_primary_gain_x; ii_stats(j).no_break_left_primary_gain_x];
    ii_results.no_break_left_primary_gain_y = [ii_results.no_break_left_primary_gain_y; ii_stats(j).no_break_left_primary_gain_y];
    ii_results.no_break_left_primary_gain_z = [ii_results.no_break_left_primary_gain_z; ii_stats(j).no_break_left_primary_gain_z];
    
    ii_results.no_break_left_final_gain_x = [ii_results.no_break_left_final_gain_x; ii_stats(j).no_break_left_final_gain_x];
    ii_results.no_break_left_final_gain_y = [ii_results.no_break_left_final_gain_y; ii_stats(j).no_break_left_final_gain_y];
    ii_results.no_break_left_final_gain_z = [ii_results.no_break_left_final_gain_z; ii_stats(j).no_break_left_final_gain_z];
    
    % SRT    
    ii_results.no_break_left_srt = [ii_results.no_break_left_srt; ii_stats(j).no_break_left_srt];
    %theta 
    %ii_results.no_break_left_primary_theta = [ii_results.no_break_left_primary_theta; ii_stats(j).no_break_left_primary_theta];
end

% Right trials
for j = 1:num_runs
    % Error
    ii_results.no_break_right_primary_err_x = [ii_results.no_break_right_primary_err_x; ii_stats(j).no_break_right_primary_err_x];
    ii_results.no_break_right_primary_err_y = [ii_results.no_break_right_primary_err_y; ii_stats(j).no_break_right_primary_err_y];
    ii_results.no_break_right_primary_err_z = [ii_results.no_break_right_primary_err_z; ii_stats(j).no_break_right_primary_err_z];

    ii_results.no_break_right_final_err_x = [ii_results.no_break_right_final_err_x; ii_stats(j).no_break_right_final_err_x];
    ii_results.no_break_right_final_err_y = [ii_results.no_break_right_final_err_y; ii_stats(j).no_break_right_final_err_y];
    ii_results.no_break_right_final_err_z = [ii_results.no_break_right_final_err_z; ii_stats(j).no_break_right_final_err_z];
    
    % Gain
    ii_results.no_break_right_primary_gain_x = [ii_results.no_break_right_primary_gain_x; ii_stats(j).no_break_right_primary_gain_x];
    ii_results.no_break_right_primary_gain_y = [ii_results.no_break_right_primary_gain_y; ii_stats(j).no_break_right_primary_gain_y];
    ii_results.no_break_right_primary_gain_z = [ii_results.no_break_right_primary_gain_z; ii_stats(j).no_break_right_primary_gain_z];
    
    ii_results.no_break_right_final_gain_x = [ii_results.no_break_right_final_gain_x; ii_stats(j).no_break_right_final_gain_x];
    ii_results.no_break_right_final_gain_y = [ii_results.no_break_right_final_gain_y; ii_stats(j).no_break_right_final_gain_y];
    ii_results.no_break_right_final_gain_z = [ii_results.no_break_right_final_gain_z; ii_stats(j).no_break_right_final_gain_z];
    
    % SRT    
    ii_results.no_break_right_srt = [ii_results.no_break_right_srt; ii_stats(j).no_break_right_srt];
     %theta 
    %ii_results.no_break_right_primary_theta = [ii_results.no_break_right_primary_theta; ii_stats(j).no_break_right_primary_theta];
end

% All trials
ii_results.mean_no_break_srt = mean(ii_results.no_break_srt);
ii_results.mean_no_break_primary_err_z = mean(ii_results.no_break_primary_err_z);
ii_results.mean_no_break_final_err_z = mean(ii_results.no_break_final_err_z);
ii_results.mean_no_break_primary_gain_z = mean(ii_results.no_break_primary_gain_z);
ii_results.mean_no_break_final_gain_z = mean(ii_results.no_break_final_gain_z);

ii_results.median_no_break_srt = median(ii_results.no_break_srt);
ii_results.median_no_break_primary_err_z = median(ii_results.no_break_primary_err_z);
ii_results.median_no_break_final_err_z = median(ii_results.no_break_final_err_z);
ii_results.median_no_break_primary_gain_z = median(ii_results.no_break_primary_gain_z);
ii_results.median_no_break_final_gain_z = median(ii_results.no_break_final_gain_z);

ii_results.std_no_break_srt = std(ii_results.no_break_srt);
ii_results.std_no_break_primary_err_z = std(ii_results.no_break_primary_err_z);
ii_results.std_no_break_final_err_z = std(ii_results.no_break_final_err_z);
ii_results.std_no_break_primary_gain_z = std(ii_results.no_break_primary_gain_z);
ii_results.std_no_break_final_gain_z = std(ii_results.no_break_final_gain_z);

% Left trials
ii_results.mean_no_break_left_srt = mean(ii_results.no_break_left_srt);
ii_results.mean_no_break_left_primary_err_z = mean(ii_results.no_break_left_primary_err_z);
ii_results.mean_no_break_left_final_err_z = mean(ii_results.no_break_left_final_err_z);
ii_results.mean_no_break_left_primary_gain_z = mean(ii_results.no_break_left_primary_gain_z);
ii_results.mean_no_break_left_final_gain_z = mean(ii_results.no_break_left_final_gain_z);

ii_results.std_no_break_left_srt = std(ii_results.no_break_left_srt);
ii_results.std_no_break_left_primary_err_z = std(ii_results.no_break_left_primary_err_z); %use this to remove trials over 2 std 


%for gg = 1:length(ii_results.no_break_left_primary_err_z);
ii_results.tossindlp =[];
     ii_results.tossindlp = find(ii_results.no_break_left_primary_err_z > 5) %3*ii_results.std_no_break_left_primary_err_z);
            newdset = removerows(ii_results.no_break_left_primary_err_z,'ind', ii_results.tossindlp);
            ii_results.no_break_left_primary_err_z_new = newdset;      
%end  


%            if ii_results.no_break_left_primary_err_z(gg) > 3*ii_results.std_no_break_left_primary_err_z;
%     ii_results.no_break_left_primary_err_z == ii_results.no_break_left_primary_err_z;
% end 
% end
ii_results.median_no_break_left_primary_err_z_new = median(ii_results.no_break_left_primary_err_z_new);
ii_results.std_no_break_left_primary_err_z_new = std(ii_results.no_break_left_primary_err_z_new);


ii_results.std_no_break_left_final_err_z = std(ii_results.no_break_left_final_err_z);
ii_results.std_no_break_left_primary_gain_z = std(ii_results.no_break_left_primary_gain_z);
ii_results.std_no_break_left_final_gain_z = std(ii_results.no_break_left_final_gain_z);

%for hh = 1:length(ii_results.no_break_left_final_err_z);
ii_results.tossindlf = [];
     ii_results.tossindlf = find(ii_results.no_break_left_final_err_z > 5) %3*ii_results.std_no_break_left_final_err_z);
            newdset = removerows(ii_results.no_break_left_final_err_z,'ind', ii_results.tossindlf);
            ii_results.no_break_left_final_err_z_new = newdset;      
%end  

ii_results.median_no_break_left_final_err_z_new = median(ii_results.no_break_left_final_err_z_new);
ii_results.std_no_break_left_final_err_z_new = std(ii_results.no_break_left_final_err_z_new);


ii_results.median_no_break_left_srt = median(ii_results.no_break_left_srt);
ii_results.median_no_break_left_primary_err_z = median(ii_results.no_break_left_primary_err_z);
ii_results.median_no_break_left_final_err_z = median(ii_results.no_break_left_final_err_z);
ii_results.median_no_break_left_primary_gain_z = median(ii_results.no_break_left_primary_gain_z);
ii_results.median_no_break_left_final_gain_z = median(ii_results.no_break_left_final_gain_z);


% Right trials
ii_results.std_no_break_right_srt = std(ii_results.no_break_right_srt);
ii_results.std_no_break_right_primary_err_z = std(ii_results.no_break_right_primary_err_z);
ii_results.std_no_break_right_final_err_z = std(ii_results.no_break_right_final_err_z);
ii_results.std_no_break_right_primary_gain_z = std(ii_results.no_break_right_primary_gain_z);
ii_results.std_no_break_right_final_gain_z = std(ii_results.no_break_right_final_gain_z);

%for rr = 1:length(ii_results.no_break_right_primary_err_z);
ii_results.tossindrp = [];
     ii_results.tossindrp = find(ii_results.no_break_right_primary_err_z > 5) %3*ii_results.std_no_break_right_primary_err_z);
            newdsetrp = removerows(ii_results.no_break_right_primary_err_z,'ind', ii_results.tossindrp);
            ii_results.no_break_right_primary_err_z_new = newdsetrp;      
%end 
%for rr = 1:length(ii_results.no_break_right_final_err_z);
ii_results.tossindrf = [];
     ii_results.tossindrf = find(ii_results.no_break_right_final_err_z > 5) %*ii_results.std_no_break_right_final_err_z);
            newdsetrf = removerows(ii_results.no_break_right_final_err_z,'ind', ii_results.tossindrf);
            ii_results.no_break_right_final_err_z_new = newdsetrf;      
%end 



ii_results.mean_no_break_right_srt = mean(ii_results.no_break_right_srt);
ii_results.mean_no_break_right_primary_err_z = mean(ii_results.no_break_right_primary_err_z);
ii_results.mean_no_break_right_final_err_z = mean(ii_results.no_break_right_final_err_z);
ii_results.mean_no_break_right_primary_gain_z = mean(ii_results.no_break_right_primary_gain_z);
ii_results.mean_no_break_right_final_gain_z = mean(ii_results.no_break_right_final_gain_z);

ii_results.median_no_break_right_srt = median(ii_results.no_break_right_srt);
ii_results.median_no_break_right_primary_err_z = median(ii_results.no_break_right_primary_err_z);
ii_results.median_no_break_right_final_err_z = median(ii_results.no_break_right_final_err_z);
ii_results.median_no_break_right_primary_gain_z = median(ii_results.no_break_right_primary_gain_z);
ii_results.median_no_break_right_final_gain_z = median(ii_results.no_break_right_final_gain_z);

ii_results.median_no_break_right_primary_err_z_new = median(ii_results.no_break_right_primary_err_z_new);
ii_results.std_no_break_right_primary_err_z_new = std(ii_results.no_break_right_primary_err_z_new);

ii_results.median_no_break_right_final_err_z_new = median(ii_results.no_break_right_final_err_z_new);
ii_results.std_no_break_right_final_err_z_new = std(ii_results.no_break_right_final_err_z_new);
%%%%%%%%%%%%%%%%%%%%%%%%
% FIXATION BREAKS ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:num_runs
    ii_stats(i).only_break_inds = find(ii_stats(i).only_break_mat==1); % Get good trials
    
    % Error
    ii_stats(i).only_break_primary_err_x = ii_stats(i).primary_err_x(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_primary_err_y = ii_stats(i).primary_err_y(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).only_break_inds);
    
    ii_stats(i).only_break_final_err_x = ii_stats(i).final_err_x(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_final_err_y = ii_stats(i).final_err_y(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_final_err_z = ii_stats(i).final_err_z(ii_stats(i).only_break_inds);
    
    % Gain
    ii_stats(i).only_break_primary_gain_x = ii_stats(i).primary_gain_x(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_primary_gain_y = ii_stats(i).primary_gain_y(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).only_break_inds);
    
    ii_stats(i).only_break_final_gain_x = ii_stats(i).final_gain_x(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_final_gain_y = ii_stats(i).final_gain_y(ii_stats(i).only_break_inds);
    ii_stats(i).only_break_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).only_break_inds);
    
    % SRT
    ii_stats(i).only_break_srt = ii_stats(i).srt(ii_stats(i).only_break_inds);
end

ii_results.only_break_primary_err_x = [];
ii_results.only_break_primary_err_y = [];
ii_results.only_break_primary_err_z = [];

ii_results.only_break_final_err_x = [];
ii_results.only_break_final_err_y = [];
ii_results.only_break_final_err_z = [];

ii_results.only_break_primary_gain_x = [];
ii_results.only_break_primary_gain_y = [];
ii_results.only_break_primary_gain_z = [];

ii_results.only_break_final_gain_x = [];
ii_results.only_break_final_gain_y = [];
ii_results.only_break_final_gain_z = [];

ii_results.only_break_srt = [];

for j = 1:num_runs
    % Error
    ii_results.only_break_primary_err_x = [ii_results.only_break_primary_err_x; ii_stats(j).only_break_primary_err_x];
    ii_results.only_break_primary_err_y = [ii_results.only_break_primary_err_y; ii_stats(j).only_break_primary_err_y];
    ii_results.only_break_primary_err_z = [ii_results.only_break_primary_err_z; ii_stats(j).only_break_primary_err_z];

    ii_results.only_break_final_err_x = [ii_results.only_break_final_err_x; ii_stats(j).only_break_final_err_x];
    ii_results.only_break_final_err_y = [ii_results.only_break_final_err_y; ii_stats(j).only_break_final_err_y];
    ii_results.only_break_final_err_z = [ii_results.only_break_final_err_z; ii_stats(j).only_break_final_err_z];
    
    % Gain
    ii_results.only_break_primary_gain_x = [ii_results.only_break_primary_gain_x; ii_stats(j).only_break_primary_gain_x];
    ii_results.only_break_primary_gain_y = [ii_results.only_break_primary_gain_y; ii_stats(j).only_break_primary_gain_y];
    ii_results.only_break_primary_gain_z = [ii_results.only_break_primary_gain_z; ii_stats(j).only_break_primary_gain_z];
    
    ii_results.only_break_final_gain_x = [ii_results.only_break_final_gain_x; ii_stats(j).only_break_final_gain_x];
    ii_results.only_break_final_gain_y = [ii_results.only_break_final_gain_y; ii_stats(j).only_break_final_gain_y];
    ii_results.only_break_final_gain_z = [ii_results.only_break_final_gain_z; ii_stats(j).only_break_final_gain_z];
    
    % SRT    
    ii_results.only_break_srt = [ii_results.only_break_srt; ii_stats(j).only_break_srt];
end

ii_results.mean_only_break_srt = mean(ii_results.only_break_srt);
ii_results.mean_only_break_primary_err_z = mean(ii_results.only_break_primary_err_z);
ii_results.mean_only_break_final_err_z = mean(ii_results.only_break_final_err_z);
ii_results.mean_only_break_primary_gain_z = mean(ii_results.only_break_primary_gain_z);
ii_results.mean_only_break_final_gain_z = mean(ii_results.only_break_final_gain_z);

ii_results.median_only_break_srt = median(ii_results.only_break_srt);
ii_results.median_only_break_primary_err_z = median(ii_results.only_break_primary_err_z);
ii_results.median_only_break_final_err_z = median(ii_results.only_break_final_err_z);
ii_results.median_only_break_primary_gain_z = median(ii_results.only_break_primary_gain_z);
ii_results.median_only_break_final_gain_z = median(ii_results.only_break_final_gain_z);

%%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT INTO HEMIFIELDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Create function to do this


%%%%%%%%%%%%%%%%%%%%%%%%
% SPLIT INTO QUADRANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Create function to do this


%%%%%%%%%%%%%%%%%%%%%%
% SHIT I WANT TO SEE %
%%%%%%%%%%%%%%%%%%%%%%

% PIE CHART FOR TRIAL COMPLIANCE
pie_nogo = ii_results.num_trials - ii_results.num_nogo;
pie_bfix = ii_results.num_bfix;
pie_else = ii_results.num_trials - pie_bfix - pie_nogo;
pie_labels = {'nogo trials', 'fixation breaks', 'good trials'};

% MAKE SUBPLOT
if cond == 'hi';
figure('Name','Trial Compliance, HIGH','NumberTitle','off')
else 
    figure('Name','Trial Compliance, LOW','NumberTitle', 'off')
end
%figure('Name','Trial compliance','NumberTitle','off')
% pie([pie_nogo, pie_bfix, pie_else]);
% legend(pie_labels);

pie([(ii_results.num_trials-pie_else),pie_else]);
% disp(ii_results.num_trials);
% disp(pie_else./ii_results.num_trials);
% disp(pie_nogo./ii_results.num_trials);
% disp(pie_bfix./ii_results.num_trials);

% DV DISTRIBUTIONS ALL
if cond == 'hi';
figure('Name','All trials, HIGH','NumberTitle','off')
else 
    figure('Name','All trials, LOW','NumberTitle', 'off')
end
%figure('Name','All Trials','NumberTitle','off')
subplot(3,2,1);
histogram(ii_results.all_left_primary_err_z,50); % 'BinWidth', (how wide bin should be)
title(['Left error: ', num2str(ii_results.median_all_left_primary_err_z), ', STD: ', num2str(ii_results.std_all_left_primary_err_z)],'FontSize',14)

subplot(3,2,2);
histogram(ii_results.all_right_primary_err_z,50);
title(['Right error: ', num2str(ii_results.median_all_right_primary_err_z), ', STD: ', num2str(ii_results.std_all_right_primary_err_z)],'FontSize',14)

subplot(3,2,3)
histogram(ii_results.all_left_primary_gain_z,50);
title(['Left gain: ', num2str(ii_results.median_all_left_primary_gain_z), ', STD: ', num2str(ii_results.std_all_left_primary_gain_z)],'FontSize',14)

subplot(3,2,4)
histogram(ii_results.all_right_primary_gain_z,50);
title(['Right gain: ', num2str(ii_results.median_all_right_primary_gain_z), ', STD: ', num2str(ii_results.std_all_right_primary_gain_z)],'FontSize',14)

subplot(3,2,5)
histogram(ii_results.all_left_srt,50);
title(['Left SRT: ', num2str(ii_results.median_all_left_srt), ', STD: ', num2str(ii_results.std_all_left_srt)],'FontSize',14)

subplot(3,2,6)
histogram(ii_results.all_right_srt,50);
title(['Right SRT: ', num2str(ii_results.median_all_right_srt), ', STD: ', num2str(ii_results.std_all_right_srt)],'FontSize',14)


% new
if cond == 'hi';
figure('Name','No break trials PRIMARY, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials PRIMARY, LOW','NumberTitle', 'off')
end
%figure('Name','No break trials PRIMARY','NumberTitle','off')
subplot(1,2,1);
histogram(ii_results.no_break_left_primary_err_z,50);
title(['Left error: ', num2str(ii_results.median_no_break_left_primary_err_z), ', STD: ', num2str(ii_results.std_no_break_left_primary_err_z)],'FontSize',14)

hold on;
subplot(1,2,2)
histogram(ii_results.no_break_left_primary_err_z_new ,50);
title(['Left error >3STD excluded: ', num2str(ii_results.median_no_break_left_primary_err_z_new), ', STD: ', num2str(ii_results.std_no_break_left_primary_err_z_new)],'FontSize',14)

if cond == 'hi';
figure('Name','No break trials PRIMARY, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials PRIMARY, LOW','NumberTitle', 'off')
end
%figure('Name','No break trials PRIMARY','NumberTitle','off')
subplot(1,2,1);
histogram(ii_results.no_break_right_primary_err_z,50);
title(['right error: ', num2str(ii_results.median_no_break_right_primary_err_z), ', STD: ', num2str(ii_results.std_no_break_right_primary_err_z)],'FontSize',14)

hold on;
subplot(1,2,2)
histogram(ii_results.no_break_right_primary_err_z_new ,50);
title(['Right error >3STD excluded: ', num2str(ii_results.median_no_break_right_primary_err_z_new), ', STD: ', num2str(ii_results.std_no_break_right_primary_err_z_new)],'FontSize',14)

% new final
if cond == 'hi';
figure('Name','No break trials FINAL, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials Final, LOW','NumberTitle', 'off')
end
%figure('Name','No break trials FINAL','NumberTitle','off')
subplot(1,2,1);
histogram(ii_results.no_break_left_final_err_z,50);
title(['Left error: ', num2str(ii_results.median_no_break_left_final_err_z), ', STD: ', num2str(ii_results.std_no_break_left_final_err_z)],'FontSize',14)

hold on;
subplot(1,2,2)
histogram(ii_results.no_break_left_final_err_z_new ,50);
title(['Left error >3STD excluded: ', num2str(ii_results.median_no_break_left_final_err_z_new), ', STD: ', num2str(ii_results.std_no_break_left_final_err_z_new)],'FontSize',14)

if cond == 'hi';
figure('Name','No break trials FINAL, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials Final, LOW','NumberTitle', 'off')
end
%figure('Name','No break trials FINAL','NumberTitle','off')
subplot(1,2,1);
histogram(ii_results.no_break_right_final_err_z,50);
title(['Right error: ', num2str(ii_results.median_no_break_right_final_err_z), ', STD: ', num2str(ii_results.std_no_break_right_final_err_z)],'FontSize',14)

hold on;
subplot(1,2,2)
histogram(ii_results.no_break_right_final_err_z_new ,50);
title(['Right error >3STD excluded:', num2str(ii_results.median_no_break_right_final_err_z_new), ', STD: ', num2str(ii_results.std_no_break_right_final_err_z_new)],'FontSize',14)

% 
% DV DISTRIBUTIONS NO BREAK PRIMARY
if cond == 'hi';
figure('Name','No break trials PRIMARY, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials PRIMARY, LOW','NumberTitle', 'off')
end
%figure('Name','No break trials PRIMARY','NumberTitle','off')
subplot(3,2,1);
histogram(ii_results.no_break_left_primary_err_z,50);
% hold on;
% histogram(ii_results.no_break_left_primary_err_z_new ,50);
title(['Left error: ', num2str(ii_results.median_no_break_left_primary_err_z), ', STD: ', num2str(ii_results.std_no_break_left_primary_err_z)],'FontSize',14)

subplot(3,2,2);
histogram(ii_results.no_break_right_primary_err_z,50);
% hold on;
% histogram(ii_results.no_break_right_primary_err_z_new,50)
title(['Right error: ', num2str(ii_results.median_no_break_right_primary_err_z), ', STD: ', num2str(ii_results.std_no_break_right_primary_err_z)],'FontSize',14)

subplot(3,2,3)
histogram(ii_results.no_break_left_primary_gain_z,50);
title(['Left gain: ', num2str(ii_results.median_no_break_left_primary_gain_z), ', STD: ', num2str(ii_results.std_no_break_left_primary_gain_z)],'FontSize',14)

subplot(3,2,4)
histogram(ii_results.no_break_right_primary_gain_z,50);
title(['Right gain: ', num2str(ii_results.median_no_break_right_primary_gain_z), ', STD: ', num2str(ii_results.std_no_break_right_primary_gain_z)],'FontSize',14)

subplot(3,2,5)
histogram(ii_results.no_break_left_srt,50);
title(['Left SRT: ', num2str(ii_results.median_no_break_left_srt), ', STD: ', num2str(ii_results.std_no_break_left_srt)],'FontSize',14)

subplot(3,2,6)
histogram(ii_results.no_break_right_srt,50);
title(['Right SRT: ', num2str(ii_results.median_no_break_right_srt), ', STD: ', num2str(ii_results.std_no_break_right_srt)],'FontSize',14)

% DV DISTRIBUTIONS NO BREAK FINAL
% DV DISTRIBUTIONS NO BREAK
if cond == 'hi';
figure('Name','No break trials FINAL, HIGH','NumberTitle','off')
else 
    figure('Name','No break trials Final, LOW','NumberTitle', 'off')
end
subplot(3,2,1);
histogram(ii_results.no_break_left_final_err_z,50);
title(['Left error: ', num2str(ii_results.median_no_break_left_final_err_z), ', STD: ', num2str(ii_results.std_no_break_left_final_err_z)],'FontSize',14)

subplot(3,2,2);
histogram(ii_results.no_break_right_final_err_z,50);
title(['Right error: ', num2str(ii_results.median_no_break_right_final_err_z), ', STD: ', num2str(ii_results.std_no_break_right_final_err_z)],'FontSize',14)

subplot(3,2,3)
histogram(ii_results.no_break_left_final_gain_z,50);
title(['Left gain: ', num2str(ii_results.median_no_break_left_final_gain_z), ', STD: ', num2str(ii_results.std_no_break_left_final_gain_z)],'FontSize',14)

subplot(3,2,4)
histogram(ii_results.no_break_right_final_gain_z,50);
title(['Right gain: ', num2str(ii_results.median_no_break_right_final_gain_z), ', STD: ', num2str(ii_results.std_no_break_right_final_gain_z)],'FontSize',14)

subplot(3,2,5)
histogram(ii_results.no_break_left_srt,50);
title(['Left SRT: ', num2str(ii_results.median_no_break_left_srt), ', STD: ', num2str(ii_results.std_no_break_left_srt)],'FontSize',14)

subplot(3,2,6)
histogram(ii_results.no_break_right_srt,50);
title(['Right SRT: ', num2str(ii_results.median_no_break_right_srt), ', STD: ', num2str(ii_results.std_no_break_right_srt)],'FontSize',14)

% ROSE PLOTS
% Create bin matrices
% Segment data into bins
% Get parameter estimates
% Plot by bin (polar(bin#,median))
% bin 1: 3.14:2.355 (2.7475)
% bin 2: 2.355:1.57 (1.9625)
% bin 3: 1.57:0.785 (1.1775)
% bin 4: 0.785:0 (0.3925)
% bin 5: 0:-0.785
% bin 6: -0.785:-1.57
% bin 7: -1.57:-2.355
% bin 8: -2.355:-3.14

for i = 1:num_runs
    
    ii_stats(i).b1_inds = find(ii_stats(i).corrective_rho > 2.355 & ii_stats(i).corrective_rho <= 3.14);
    ii_stats(i).b2_inds = find(ii_stats(i).corrective_rho > 1.57 & ii_stats(i).corrective_rho <= 2.355);
    ii_stats(i).b3_inds = find(ii_stats(i).corrective_rho > 0.785 & ii_stats(i).corrective_rho <= 1.57);
    ii_stats(i).b4_inds = find(ii_stats(i).corrective_rho > 0 & ii_stats(i).corrective_rho <= 0.785);
    ii_stats(i).b5_inds = find(ii_stats(i).corrective_rho > -0.785 & ii_stats(i).corrective_rho <= 0);
    ii_stats(i).b6_inds = find(ii_stats(i).corrective_rho > -1.57 & ii_stats(i).corrective_rho <= -0.785);
    ii_stats(i).b7_inds = find(ii_stats(i).corrective_rho > -2.355 & ii_stats(i).corrective_rho <= -1.57);
    ii_stats(i).b8_inds = find(ii_stats(i).corrective_rho > -3.14 & ii_stats(i).corrective_rho <= -2.355);

    ii_stats(i).b1_mat = zeros(trialnum,1);
    ii_stats(i).b2_mat = zeros(trialnum,1);
    ii_stats(i).b3_mat = zeros(trialnum,1);
    ii_stats(i).b4_mat = zeros(trialnum,1);
    ii_stats(i).b5_mat = zeros(trialnum,1);
    ii_stats(i).b6_mat = zeros(trialnum,1);
    ii_stats(i).b7_mat = zeros(trialnum,1);
    ii_stats(i).b8_mat = zeros(trialnum,1);
    
    ii_stats(i).b1_mat(ii_stats(i).b1_inds) = 1;
    ii_stats(i).b2_mat(ii_stats(i).b2_inds) = 1;
    ii_stats(i).b3_mat(ii_stats(i).b3_inds) = 1;
    ii_stats(i).b4_mat(ii_stats(i).b4_inds) = 1;
    ii_stats(i).b5_mat(ii_stats(i).b5_inds) = 1;
    ii_stats(i).b6_mat(ii_stats(i).b6_inds) = 1;
    ii_stats(i).b7_mat(ii_stats(i).b7_inds) = 1;
    ii_stats(i).b8_mat(ii_stats(i).b8_inds) = 1;
    
    ii_stats(i).no_break_b1_mat = ii_stats(i).no_break_mat .* ii_stats(i).b1_mat;
    ii_stats(i).no_break_b2_mat = ii_stats(i).no_break_mat .* ii_stats(i).b2_mat;
    ii_stats(i).no_break_b3_mat = ii_stats(i).no_break_mat .* ii_stats(i).b3_mat;
    ii_stats(i).no_break_b4_mat = ii_stats(i).no_break_mat .* ii_stats(i).b4_mat;
    ii_stats(i).no_break_b5_mat = ii_stats(i).no_break_mat .* ii_stats(i).b5_mat;
    ii_stats(i).no_break_b6_mat = ii_stats(i).no_break_mat .* ii_stats(i).b6_mat;
    ii_stats(i).no_break_b7_mat = ii_stats(i).no_break_mat .* ii_stats(i).b7_mat;
    ii_stats(i).no_break_b8_mat = ii_stats(i).no_break_mat .* ii_stats(i).b8_mat;
    
end

for i = 1:num_runs
    
    ii_stats(i).no_break_b1_inds = find(ii_stats(i).no_break_b1_mat==1);
    ii_stats(i).no_break_b2_inds = find(ii_stats(i).no_break_b2_mat==1);
    ii_stats(i).no_break_b3_inds = find(ii_stats(i).no_break_b3_mat==1);
    ii_stats(i).no_break_b4_inds = find(ii_stats(i).no_break_b4_mat==1);
    ii_stats(i).no_break_b5_inds = find(ii_stats(i).no_break_b5_mat==1);
    ii_stats(i).no_break_b6_inds = find(ii_stats(i).no_break_b6_mat==1);
    ii_stats(i).no_break_b7_inds = find(ii_stats(i).no_break_b7_mat==1);
    ii_stats(i).no_break_b8_inds = find(ii_stats(i).no_break_b8_mat==1);
    
    % PRIMARY ERROR
    ii_stats(i).no_break_b1_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b1_inds);
    ii_stats(i).no_break_b2_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b2_inds);
    ii_stats(i).no_break_b3_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b3_inds);
    ii_stats(i).no_break_b4_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b4_inds);
    ii_stats(i).no_break_b5_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b5_inds);
    ii_stats(i).no_break_b6_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b6_inds);
    ii_stats(i).no_break_b7_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b7_inds);
    ii_stats(i).no_break_b8_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).no_break_b8_inds);
    
    % PRIMARY GAIN
    ii_stats(i).no_break_b1_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b1_inds);
    ii_stats(i).no_break_b2_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b2_inds);
    ii_stats(i).no_break_b3_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b3_inds);
    ii_stats(i).no_break_b4_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b4_inds);
    ii_stats(i).no_break_b5_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b5_inds);
    ii_stats(i).no_break_b6_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b6_inds);
    ii_stats(i).no_break_b7_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b7_inds);
    ii_stats(i).no_break_b8_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).no_break_b8_inds);
    
    % SRT
    ii_stats(i).no_break_b1_srt = ii_stats(i).srt(ii_stats(i).no_break_b1_inds);
    ii_stats(i).no_break_b2_srt = ii_stats(i).srt(ii_stats(i).no_break_b2_inds);
    ii_stats(i).no_break_b3_srt = ii_stats(i).srt(ii_stats(i).no_break_b3_inds);
    ii_stats(i).no_break_b4_srt = ii_stats(i).srt(ii_stats(i).no_break_b4_inds);
    ii_stats(i).no_break_b5_srt = ii_stats(i).srt(ii_stats(i).no_break_b5_inds);
    ii_stats(i).no_break_b6_srt = ii_stats(i).srt(ii_stats(i).no_break_b6_inds);
    ii_stats(i).no_break_b7_srt = ii_stats(i).srt(ii_stats(i).no_break_b7_inds);
    ii_stats(i).no_break_b8_srt = ii_stats(i).srt(ii_stats(i).no_break_b8_inds);

    % FINAL ERROR
    ii_stats(i).no_break_b1_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b1_inds);
    ii_stats(i).no_break_b2_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b2_inds);
    ii_stats(i).no_break_b3_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b3_inds);
    ii_stats(i).no_break_b4_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b4_inds);
    ii_stats(i).no_break_b5_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b5_inds);
    ii_stats(i).no_break_b6_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b6_inds);
    ii_stats(i).no_break_b7_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b7_inds);
    ii_stats(i).no_break_b8_final_err_z = ii_stats(i).final_err_z(ii_stats(i).no_break_b8_inds);
    
    % FINAL GAIN
    ii_stats(i).no_break_b1_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b1_inds);
    ii_stats(i).no_break_b2_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b2_inds);
    ii_stats(i).no_break_b3_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b3_inds);
    ii_stats(i).no_break_b4_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b4_inds);
    ii_stats(i).no_break_b5_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b5_inds);
    ii_stats(i).no_break_b6_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b6_inds);
    ii_stats(i).no_break_b7_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b7_inds);
    ii_stats(i).no_break_b8_final_gain_z = ii_stats(i).final_gain_z(ii_stats(i).no_break_b8_inds);
 
%     % PRIMARY ERROR
%     ii_stats(i).b1_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b1_inds);
%     ii_stats(i).b2_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b2_inds);
%     ii_stats(i).b3_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b3_inds);
%     ii_stats(i).b4_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b4_inds);
%     ii_stats(i).b5_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b5_inds);
%     ii_stats(i).b6_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b6_inds);
%     ii_stats(i).b7_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b7_inds);
%     ii_stats(i).b8_primary_err_z = ii_stats(i).primary_err_z(ii_stats(i).b8_inds);
%     
%     % PRIMARY GAIN
%     ii_stats(i).b1_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b1_inds);
%     ii_stats(i).b2_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b2_inds);
%     ii_stats(i).b3_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b3_inds);
%     ii_stats(i).b4_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b4_inds);
%     ii_stats(i).b5_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b5_inds);
%     ii_stats(i).b6_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b6_inds);
%     ii_stats(i).b7_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b7_inds);
%     ii_stats(i).b8_primary_gain_z = ii_stats(i).primary_gain_z(ii_stats(i).b8_inds);
%     
%     % SRT
%     ii_stats(i).b1_srt = ii_stats(i).srt(ii_stats(i).b1_inds);
%     ii_stats(i).b2_srt = ii_stats(i).srt(ii_stats(i).b2_inds);
%     ii_stats(i).b3_srt = ii_stats(i).srt(ii_stats(i).b3_inds);
%     ii_stats(i).b4_srt = ii_stats(i).srt(ii_stats(i).b4_inds);
%     ii_stats(i).b5_srt = ii_stats(i).srt(ii_stats(i).b5_inds);
%     ii_stats(i).b6_srt = ii_stats(i).srt(ii_stats(i).b6_inds);
%     ii_stats(i).b7_srt = ii_stats(i).srt(ii_stats(i).b7_inds);
%     ii_stats(i).b8_srt = ii_stats(i).srt(ii_stats(i).b8_inds);
    
end

ii_results.b1_all_primary_err_z = [];
ii_results.b2_all_primary_err_z = [];
ii_results.b3_all_primary_err_z = [];
ii_results.b4_all_primary_err_z = [];
ii_results.b5_all_primary_err_z = [];
ii_results.b6_all_primary_err_z = [];
ii_results.b7_all_primary_err_z = [];
ii_results.b8_all_primary_err_z = [];

ii_results.b1_all_primary_gain_z = [];
ii_results.b2_all_primary_gain_z = [];
ii_results.b3_all_primary_gain_z = [];
ii_results.b4_all_primary_gain_z = [];
ii_results.b5_all_primary_gain_z = [];
ii_results.b6_all_primary_gain_z = [];
ii_results.b7_all_primary_gain_z = [];
ii_results.b8_all_primary_gain_z = [];

ii_results.b1_all_srt = [];
ii_results.b2_all_srt = [];
ii_results.b3_all_srt = [];
ii_results.b4_all_srt = [];
ii_results.b5_all_srt = [];
ii_results.b6_all_srt = [];
ii_results.b7_all_srt = [];
ii_results.b8_all_srt = [];

ii_results.b1_no_break_primary_err_z = [];
ii_results.b2_no_break_primary_err_z = [];
ii_results.b3_no_break_primary_err_z = [];
ii_results.b4_no_break_primary_err_z = [];
ii_results.b5_no_break_primary_err_z = [];
ii_results.b6_no_break_primary_err_z = [];
ii_results.b7_no_break_primary_err_z = [];
ii_results.b8_no_break_primary_err_z = [];

ii_results.b1_no_break_primary_gain_z = [];
ii_results.b2_no_break_primary_gain_z = [];
ii_results.b3_no_break_primary_gain_z = [];
ii_results.b4_no_break_primary_gain_z = [];
ii_results.b5_no_break_primary_gain_z = [];
ii_results.b6_no_break_primary_gain_z = [];
ii_results.b7_no_break_primary_gain_z = [];
ii_results.b8_no_break_primary_gain_z = [];

ii_results.b1_no_break_srt = [];
ii_results.b2_no_break_srt = [];
ii_results.b3_no_break_srt = [];
ii_results.b4_no_break_srt = [];
ii_results.b5_no_break_srt = [];
ii_results.b6_no_break_srt = [];
ii_results.b7_no_break_srt = [];
ii_results.b8_no_break_srt = [];

ii_results.b1_no_break_final_err_z = [];
ii_results.b2_no_break_final_err_z = [];
ii_results.b3_no_break_final_err_z = [];
ii_results.b4_no_break_final_err_z = [];
ii_results.b5_no_break_final_err_z = [];
ii_results.b6_no_break_final_err_z = [];
ii_results.b7_no_break_final_err_z = [];
ii_results.b8_no_break_final_err_z = [];

ii_results.b1_no_break_final_gain_z = [];
ii_results.b2_no_break_final_gain_z = [];
ii_results.b3_no_break_final_gain_z = [];
ii_results.b4_no_break_final_gain_z = [];
ii_results.b5_no_break_final_gain_z = [];
ii_results.b6_no_break_final_gain_z = [];
ii_results.b7_no_break_final_gain_z = [];
ii_results.b8_no_break_final_gain_z = [];

for j = 1:num_runs
    % Error
%     ii_results.b1_all_primary_err_z = [ii_results.b1_all_primary_err_z; ii_stats(j).b1_primary_err_z];
%     ii_results.b2_all_primary_err_z = [ii_results.b2_all_primary_err_z; ii_stats(j).b2_primary_err_z];
%     ii_results.b3_all_primary_err_z = [ii_results.b3_all_primary_err_z; ii_stats(j).b3_primary_err_z];
%     ii_results.b4_all_primary_err_z = [ii_results.b4_all_primary_err_z; ii_stats(j).b4_primary_err_z];
%     ii_results.b5_all_primary_err_z = [ii_results.b5_all_primary_err_z; ii_stats(j).b5_primary_err_z];
%     ii_results.b6_all_primary_err_z = [ii_results.b6_all_primary_err_z; ii_stats(j).b6_primary_err_z];
%     ii_results.b7_all_primary_err_z = [ii_results.b7_all_primary_err_z; ii_stats(j).b7_primary_err_z];
%     ii_results.b8_all_primary_err_z = [ii_results.b8_all_primary_err_z; ii_stats(j).b8_primary_err_z];
%     
%     % Gain
%     ii_results.b1_all_primary_gain_z = [ii_results.b1_all_primary_gain_z; ii_stats(j).b1_primary_gain_z];
%     ii_results.b2_all_primary_gain_z = [ii_results.b2_all_primary_gain_z; ii_stats(j).b2_primary_gain_z];
%     ii_results.b3_all_primary_gain_z = [ii_results.b3_all_primary_gain_z; ii_stats(j).b3_primary_gain_z];
%     ii_results.b4_all_primary_gain_z = [ii_results.b4_all_primary_gain_z; ii_stats(j).b4_primary_gain_z];
%     ii_results.b5_all_primary_gain_z = [ii_results.b5_all_primary_gain_z; ii_stats(j).b5_primary_gain_z];
%     ii_results.b6_all_primary_gain_z = [ii_results.b6_all_primary_gain_z; ii_stats(j).b6_primary_gain_z];
%     ii_results.b7_all_primary_gain_z = [ii_results.b7_all_primary_gain_z; ii_stats(j).b7_primary_gain_z];
%     ii_results.b8_all_primary_gain_z = [ii_results.b8_all_primary_gain_z; ii_stats(j).b8_primary_gain_z];
%     
%     % SRT    
%     ii_results.b1_all_srt = [ii_results.b1_all_srt; ii_stats(j).b1_srt];
%     ii_results.b2_all_srt = [ii_results.b2_all_srt; ii_stats(j).b2_srt];
%     ii_results.b3_all_srt = [ii_results.b3_all_srt; ii_stats(j).b3_srt];
%     ii_results.b4_all_srt = [ii_results.b4_all_srt; ii_stats(j).b4_srt];
%     ii_results.b5_all_srt = [ii_results.b5_all_srt; ii_stats(j).b5_srt];
%     ii_results.b6_all_srt = [ii_results.b6_all_srt; ii_stats(j).b6_srt];
%     ii_results.b7_all_srt = [ii_results.b7_all_srt; ii_stats(j).b7_srt];
%     ii_results.b8_all_srt = [ii_results.b8_all_srt; ii_stats(j).b8_srt];

    ii_results.b1_no_break_primary_err_z = [ii_results.b1_no_break_primary_err_z; ii_stats(j).no_break_b1_primary_err_z];
    ii_results.b2_no_break_primary_err_z = [ii_results.b2_no_break_primary_err_z; ii_stats(j).no_break_b2_primary_err_z];
    ii_results.b3_no_break_primary_err_z = [ii_results.b3_no_break_primary_err_z; ii_stats(j).no_break_b3_primary_err_z];
    ii_results.b4_no_break_primary_err_z = [ii_results.b4_no_break_primary_err_z; ii_stats(j).no_break_b4_primary_err_z];
    ii_results.b5_no_break_primary_err_z = [ii_results.b5_no_break_primary_err_z; ii_stats(j).no_break_b5_primary_err_z];
    ii_results.b6_no_break_primary_err_z = [ii_results.b6_no_break_primary_err_z; ii_stats(j).no_break_b6_primary_err_z];
    ii_results.b7_no_break_primary_err_z = [ii_results.b7_no_break_primary_err_z; ii_stats(j).no_break_b7_primary_err_z];
    ii_results.b8_no_break_primary_err_z = [ii_results.b8_no_break_primary_err_z; ii_stats(j).no_break_b8_primary_err_z];
    
    % Gain
    ii_results.b1_no_break_primary_gain_z = [ii_results.b1_no_break_primary_gain_z; ii_stats(j).no_break_b1_primary_gain_z];
    ii_results.b2_no_break_primary_gain_z = [ii_results.b2_no_break_primary_gain_z; ii_stats(j).no_break_b2_primary_gain_z];
    ii_results.b3_no_break_primary_gain_z = [ii_results.b3_no_break_primary_gain_z; ii_stats(j).no_break_b3_primary_gain_z];
    ii_results.b4_no_break_primary_gain_z = [ii_results.b4_no_break_primary_gain_z; ii_stats(j).no_break_b4_primary_gain_z];
    ii_results.b5_no_break_primary_gain_z = [ii_results.b5_no_break_primary_gain_z; ii_stats(j).no_break_b5_primary_gain_z];
    ii_results.b6_no_break_primary_gain_z = [ii_results.b6_no_break_primary_gain_z; ii_stats(j).no_break_b6_primary_gain_z];
    ii_results.b7_no_break_primary_gain_z = [ii_results.b7_no_break_primary_gain_z; ii_stats(j).no_break_b7_primary_gain_z];
    ii_results.b8_no_break_primary_gain_z = [ii_results.b8_no_break_primary_gain_z; ii_stats(j).no_break_b8_primary_gain_z];
    
    % SRT    
    ii_results.b1_no_break_srt = [ii_results.b1_no_break_srt; ii_stats(j).no_break_b1_srt];
    ii_results.b2_no_break_srt = [ii_results.b2_no_break_srt; ii_stats(j).no_break_b2_srt];
    ii_results.b3_no_break_srt = [ii_results.b3_no_break_srt; ii_stats(j).no_break_b3_srt];
    ii_results.b4_no_break_srt = [ii_results.b4_no_break_srt; ii_stats(j).no_break_b4_srt];
    ii_results.b5_no_break_srt = [ii_results.b5_no_break_srt; ii_stats(j).no_break_b5_srt];
    ii_results.b6_no_break_srt = [ii_results.b6_no_break_srt; ii_stats(j).no_break_b6_srt];
    ii_results.b7_no_break_srt = [ii_results.b7_no_break_srt; ii_stats(j).no_break_b7_srt];
    ii_results.b8_no_break_srt = [ii_results.b8_no_break_srt; ii_stats(j).no_break_b8_srt];
    
    % Final Error
    ii_results.b1_no_break_final_err_z = [ii_results.b1_no_break_final_err_z; ii_stats(j).no_break_b1_final_err_z];
    ii_results.b2_no_break_final_err_z = [ii_results.b2_no_break_final_err_z; ii_stats(j).no_break_b2_final_err_z];
    ii_results.b3_no_break_final_err_z = [ii_results.b3_no_break_final_err_z; ii_stats(j).no_break_b3_final_err_z];
    ii_results.b4_no_break_final_err_z = [ii_results.b4_no_break_final_err_z; ii_stats(j).no_break_b4_final_err_z];
    ii_results.b5_no_break_final_err_z = [ii_results.b5_no_break_final_err_z; ii_stats(j).no_break_b5_final_err_z];
    ii_results.b6_no_break_final_err_z = [ii_results.b6_no_break_final_err_z; ii_stats(j).no_break_b6_final_err_z];
    ii_results.b7_no_break_final_err_z = [ii_results.b7_no_break_final_err_z; ii_stats(j).no_break_b7_final_err_z];
    ii_results.b8_no_break_final_err_z = [ii_results.b8_no_break_final_err_z; ii_stats(j).no_break_b8_final_err_z];
    
    % Final Gain
    ii_results.b1_no_break_final_gain_z = [ii_results.b1_no_break_final_gain_z; ii_stats(j).no_break_b1_final_gain_z];
    ii_results.b2_no_break_final_gain_z = [ii_results.b2_no_break_final_gain_z; ii_stats(j).no_break_b2_final_gain_z];
    ii_results.b3_no_break_final_gain_z = [ii_results.b3_no_break_final_gain_z; ii_stats(j).no_break_b3_final_gain_z];
    ii_results.b4_no_break_final_gain_z = [ii_results.b4_no_break_final_gain_z; ii_stats(j).no_break_b4_final_gain_z];
    ii_results.b5_no_break_final_gain_z = [ii_results.b5_no_break_final_gain_z; ii_stats(j).no_break_b5_final_gain_z];
    ii_results.b6_no_break_final_gain_z = [ii_results.b6_no_break_final_gain_z; ii_stats(j).no_break_b6_final_gain_z];
    ii_results.b7_no_break_final_gain_z = [ii_results.b7_no_break_final_gain_z; ii_stats(j).no_break_b7_final_gain_z];
    ii_results.b8_no_break_final_gain_z = [ii_results.b8_no_break_final_gain_z; ii_stats(j).no_break_b8_final_gain_z];


end

bins = [2.7475; 1.9625; 1.1775; 0.3925; -0.3925; -1.1775; -1.9625; -2.7475;2.7475];
% dt_err = [median(ii_results.b1_all_primary_err_z); median(ii_results.b2_all_primary_err_z); median(ii_results.b3_all_primary_err_z); median(ii_results.b4_all_primary_err_z); median(ii_results.b5_all_primary_err_z); median(ii_results.b6_all_primary_err_z); median(ii_results.b7_all_primary_err_z); median(ii_results.b8_all_primary_err_z)];
% dt_gain = [median(ii_results.b1_all_primary_gain_z); median(ii_results.b2_all_primary_gain_z); median(ii_results.b3_all_primary_gain_z); median(ii_results.b4_all_primary_gain_z); median(ii_results.b5_all_primary_gain_z); median(ii_results.b6_all_primary_gain_z); median(ii_results.b7_all_primary_gain_z); median(ii_results.b8_all_primary_gain_z)];
% dt_srt = [median(ii_results.b1_all_srt); median(ii_results.b2_all_srt); median(ii_results.b3_all_srt); median(ii_results.b4_all_srt); median(ii_results.b5_all_srt); median(ii_results.b6_all_srt); median(ii_results.b7_all_srt); median(ii_results.b8_all_srt)];

dt_err = [median(ii_results.b1_no_break_primary_err_z); median(ii_results.b2_no_break_primary_err_z); median(ii_results.b3_no_break_primary_err_z); median(ii_results.b4_no_break_primary_err_z); median(ii_results.b5_no_break_primary_err_z); median(ii_results.b6_no_break_primary_err_z); median(ii_results.b7_no_break_primary_err_z); median(ii_results.b8_no_break_primary_err_z);median(ii_results.b1_no_break_primary_err_z)];
dt_gain = [median(ii_results.b1_no_break_primary_gain_z); median(ii_results.b2_no_break_primary_gain_z); median(ii_results.b3_no_break_primary_gain_z); median(ii_results.b4_no_break_primary_gain_z); median(ii_results.b5_no_break_primary_gain_z); median(ii_results.b6_no_break_primary_gain_z); median(ii_results.b7_no_break_primary_gain_z); median(ii_results.b8_no_break_primary_gain_z);median(ii_results.b1_no_break_primary_gain_z)];
dt_srt = [median(ii_results.b1_no_break_srt); median(ii_results.b2_no_break_srt); median(ii_results.b3_no_break_srt); median(ii_results.b4_no_break_srt); median(ii_results.b5_no_break_srt); median(ii_results.b6_no_break_srt); median(ii_results.b7_no_break_srt); median(ii_results.b8_no_break_srt);median(ii_results.b1_no_break_srt)];

dt_err_f = [median(ii_results.b1_no_break_final_err_z); median(ii_results.b2_no_break_final_err_z); median(ii_results.b3_no_break_final_err_z); median(ii_results.b4_no_break_final_err_z); median(ii_results.b5_no_break_final_err_z); median(ii_results.b6_no_break_final_err_z); median(ii_results.b7_no_break_final_err_z); median(ii_results.b8_no_break_final_err_z);median(ii_results.b1_no_break_final_err_z)];
dt_gain_f = [median(ii_results.b1_no_break_final_gain_z); median(ii_results.b2_no_break_final_gain_z); median(ii_results.b3_no_break_final_gain_z); median(ii_results.b4_no_break_final_gain_z); median(ii_results.b5_no_break_final_gain_z); median(ii_results.b6_no_break_final_gain_z); median(ii_results.b7_no_break_final_gain_z); median(ii_results.b8_no_break_final_gain_z);median(ii_results.b1_no_break_final_gain_z)];


% figure('Name','No break trials primary error','NumberTitle','off')
% polar(bins,dt_err,'*');
figure('Name','No break trials primary error','NumberTitle','off')
p = mmpolar(bins,dt_err,'-ko','RLimit',[0 5],'TTickDelta',45);

figure('Name','No break trials primary gain','NumberTitle','off')
%polar(bins,dt_gain,'*');
p = mmpolar(bins,dt_gain,'-ko','RLimit',[0 1.2],'TTickDelta',45);

figure('Name','No break trials SRT','NumberTitle','off')
%polar(bins,dt_srt,'*');
p = mmpolar(bins,dt_srt,'-ko','RLimit',[0 500],'TTickDelta',45);

figure('Name','No break trials final error','NumberTitle','off')
p = mmpolar(bins,dt_err_f,'-ko','RLimit',[0 5],'TTickDelta',45);

figure('Name','No break trials final gain','NumberTitle','off')
%polar(bins,dt_gain,'*');
p = mmpolar(bins,dt_gain_f,'-ko','RLimit',[0 1.2],'TTickDelta',45);

if cond == 'hi'
    ii_results_hi = ii_results; 
    save('ii_results_hi')
else
ii_results_lo = ii_results; 
save('ii_results_lo')
end 
 clear ii_stats.mat
end
%putvar(bins,dt_err,dt_err_f); 04/06/17 is this absolutely necessary?

