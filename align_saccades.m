%% align all targets to common reference point
%SUMMARY: this script takes saccadic endpoint data from trials and
%normalizes them to a single theta value with rho of variable length.
%it takes the target location from ii_cfg.trialvec and operates on all X
%and Y values to transform them to the upper right quadrant of the screen.
%then, those X/Y values are converted to polar coordinates. If the theta
%component is greater/less than than 45 degrees (in radians), 45deginrad is
%subtracted/added from/to theta to acquire a difference in radians from 45.Then,
%this difference is subtracted from the true saccadic endpoint theta
%(saccadic endpoints have also been converted from cartesian to polar
%coordinates). In so doing, the theta related variance in the actual eye
%movement is preserved, while a distance is subtracted such that the eye
%movement is essentially realigned to different point in space. 
addpath(genpath('/Volumes/hyper/experiments/Grace/iEye_irresponsible-master'))

%%
load('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/ii_stats')
theta_vect=[];
rho_vect =[];
all_theta45 =[];
all_Tar_theta =[];
all_TarX =[];
all_TarY =[];
all_Tar_rho =[];
%runs = {'run02','run03','run04','run05'};
runs ={'run11','run12','run22','run23','run24','run25','run26','run27','run28'};
%runs ={'run02'};

priority_id = [];
all_newx =  [];
all_newy =  [];
all_al=[];
%runnum ={'2','3','4','5'};
runnum= {'11','12','22','23','24','25','26','27','28'};

for jj = 1:length(runs);
filename = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/TASK/%s.mat',runs{jj})
fileID = load(filename)
if  jj == 1;
    newrow = 1;
    endrow = 36;
else
newrow = ((jj-1).*36)+1;
endrow = newrow + 35;
end; 
priority_id(newrow:endrow,1) = fileID.task.conditionAndQueriedTarget(:,1);

runnumtmp = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/DATA_PROC/mr_1_%s_proc.mat',runnum{jj});
load(runnumtmp)
figure;%iEye;
ii_definetrial('XDAT',1,'XDAT',6)

for ii = 1:36 %:36 %[1 24];
tmp = [ii_cfg.trialvec XDAT]; %concatenate necessary vecs of equal length
test = tmp(:,1)==(ii) & tmp(:,2)==4;
XPri = X(test);
YPri = Y(test);
test2 = tmp(:,1)==(ii) & tmp(:,2)==5;
Xfinal = X(test2);
TarXnew = TarX(test2);
all_TarX = [all_TarX; TarXnew(1)];
TarYnew = TarY(test2);
all_TarY = [all_TarY; TarYnew(1)];
[Tar_theta Tar_rho] = cart2pol(TarXnew,TarYnew);
all_Tar_theta = [all_Tar_theta; Tar_theta(1)]; 
all_Tar_rho =[all_Tar_rho; Tar_rho(1)];

Tar_theta = Tar_theta(1);

if TarXnew(end) <0
  TarXnew = TarXnew(end) .*-1;
else
 TarXnew = TarXnew(end);
end

if TarYnew (end) < 0
   TarYnew  = TarYnew(end) .*-1;
else
   TarYnew  = TarYnew(end);
end

if XPri(end) <0
    XPriend = XPri(end) .*-1;
else
    XPriend = XPri(end);
end

if YPri(end) < 0
    YPriend = YPri(end) .*-1;
else
    YPriend = YPri(end);
end

fortyfivedeginrad = 0.785318; 
  
    if Tar_theta <  fortyfivedeginrad; 
        
         diff = fortyfivedeginrad - Tar_theta; 
         all_Tar_theta_al = Tar_theta + diff;
              all_al =[all_al; all_Tar_theta_al];
    else
       Tar_theta > fortyfivedeginrad
        diff = Tar_theta  - fortyfivedeginrad;
        all_Tar_theta_al = Tar_theta - diff;
        all_al =[all_al; all_Tar_theta_al];
    end
    

[Tar_theta Tar_rho] = cart2pol(TarXnew,TarYnew);

[Priend_theta Priend_rho] = cart2pol(XPriend,YPriend);
theta_vect =[theta_vect; Priend_theta];
rho_vect = [rho_vect; Priend_rho];

fortyfivedeginrad = 0.785318; 
  
    if Tar_theta <  fortyfivedeginrad; 
        
         diff = fortyfivedeginrad - Tar_theta; 
         theta_plot = Priend_theta + diff;
    else
       Tar_theta > fortyfivedeginrad;
        diff = Tar_theta  - fortyfivedeginrad;
        theta_plot = Priend_theta  - diff;
    end
    
all_theta45 = [all_theta45; theta_plot];

 [newx, newy]= pol2cart(theta_plot,Priend_rho);
 all_newx =  [all_newx; newx];
  all_newy =  [all_newy; newy];

if priority_id(ii) == 31; 
    newxhigh =newx;
    newyhigh =newy;
    figure(2);
plot(newxhigh, newyhigh,'bo','markersize',4)%,'markersize',4) %bo','markersize',4)
 xlim([-6 12])
 ylim([-6 12])  
hold on;
    
else 
     newxlow =newx;
    newylow =newy;
    %figure(3)
plot(newxlow,newylow,'ro','markersize',4)
 xlim([-6 12])
 ylim([-6 12]) 
hold on;
end

end 
end
figure(4);
plot(abs(all_TarX),abs(all_TarY),'k*')
figure(4)
polar(all_al, all_Tar_rho,'ko')

%%
%conferred with mr, says i can treat rho dimension as same as theta dim.
% will take a delta btwn x,y from reference point of 
load('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/ii_stats')

all_TarX =[];
all_TarY =[];
all_Tar_rho =[];
runs ={'run11','run12','run22','run23','run24','run25','run26','run27','run28'};

all_aligny = [];
all_alignx = [];
all_XPriend = []; 
all_YPriend = []; 

runnum= {'11','12','22','23','24','25','26','27','28'};

all_priority_id = [];
for jj = 1:length(runs);
    priority_id = [];
    filename = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/TASK/%s.mat',runs{jj})
    fileID = load(filename)
    priority_id = [priority_id; fileID.task.conditionAndQueriedTarget(:,1);];
    all_priority_id = [all_priority_id; priority_id]
    runnumtmp = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/DATA_PROC/mr_1_%s_proc.mat',runnum{jj});
    load(runnumtmp)
    
    ii_definetrial('XDAT',1,'XDAT',6)

for ii = 1:36 %:36 %[1 24];
tmp = [ii_cfg.trialvec XDAT]; %concatenate necessary vecs of equal length
test = tmp(:,1)==(ii) & tmp(:,2)==4;
XPri = X(test);
YPri = Y(test);
test2 = tmp(:,1)==(ii) & tmp(:,2)==5;
Xfinal = X(test2);
TarXnew = TarX(test2);
all_TarX = [all_TarX; TarXnew(1)];
TarYnew = TarY(test2);
all_TarY = [all_TarY; TarYnew(1)];



TarXnew = TarXnew(1);
TarYnew =TarYnew(1);

if TarXnew <0
  TarXnew = TarXnew .*-1;
else
 TarXnew = TarXnew;
end

if TarYnew < 0
   TarYnew  = TarYnew .*-1;
else
   TarYnew  = TarYnew;
end
Xone = XPri(1);
Yone = YPri(1);

if XPri(end) < 0
    XPriend = XPri(end) .*-1;% these are trials which occurred on the left 


else
    XPriend = XPri(end);
end
all_XPriend = [all_XPriend; Xone]; 
if YPri(end) < 0
    YPriend = YPri(end) .*-1;
else
    YPriend = YPri(end);
end

all_YPriend = [all_YPriend; Yone]; 

avgx = 5.9685; %an avg of the runs of this ppt 
avgy = 6.1394; 
  
    if TarXnew <  avgx; 
         diffx = avgx - TarXnew;          
    else
        
        TarYnew > avgy
        diffx = TarYnew  - avgy;
    end
    
       if TarYnew <  avgy; 
         diffy = avgx - TarYnew;          
    else
        
        TarYnew > avgy;
        diffy = TarYnew  - avgy;
    end
    
    
    if XPriend <  avgx; 
         alignx = XPriend + diffx;
    else
       XPriend > avgx; 
        alignx = XPriend - diffx;
    end
    
        if XPriend <  avgx; 
         alignx = XPriend + diffx;
    else
       XPriend > avgx; 
        alignx = XPriend - diffx;
    end
    
    if YPriend <  avgy; 
         aligny = YPriend + diffy;
    else
       YPriend > avgy; 
        aligny = YPriend - diffy;
    end
    
    all_alignx =  [all_alignx; alignx];
    all_aligny =  [all_aligny; aligny];
end 
end

%% 
all_left = find(all_XPriend < 0); %these indices will match aligned x&y
all_right = find(all_XPriend > 0); 
x_high = [];
y_high = [];
x_low = [];
y_low= [];


priority_id_high = find(all_priority_id == 31);
high  = all_XPriend(priority_id_high);
highlefttmp = find(high < 0);
highrighttmp = find(high > 0);
x_left_high =  all_alignx(highlefttmp);
y_left_high = all_aligny(highlefttmp);

x_right_high =  all_alignx(highrighttmp);
y_right_high = all_aligny(highrighttmp);

      
priority_id_low = find(all_priority_id == 32);
low  = all_XPriend(priority_id_low);
lowlefttmp = find(low < 0);
lowrighttmp = find(low > 0);

x_left_low =  all_alignx(lowlefttmp);
y_left_low = all_aligny(lowlefttmp);
x_right_low =  all_alignx(lowrighttmp);
y_right_low = all_aligny(lowrighttmp);

figure(2);
plot(x_left_high,y_left_high,'b*','markersize',3)%,'markersize',4) %bo','markersize',4) hold on;
hold on;
plot(x_right_high,y_right_high,'go','markersize',3)%,'markersize',4) %bo','markersize',4)
plot(x_left_low,y_left_low,'r*','markersize',3)%,'markersize',4) %bo','markersize',4)
plot(x_right_low,y_right_low,'k*','markersize',3)%
hold on;
plot(avgx, avgy, 'k+','markersize',10)

%% visualize with targ + err

load('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/ii_stats')

all_TarX =[];
all_TarY =[];
all_Tar_rho =[];
runs ={'run11','run12','run22','run23','run24','run25','run26','run27','run28'};

all_aligny = [];
all_alignx = [];
all_XPriend = []; 
all_YPriend = []; 

runnum= {'11','12','22','23','24','25','26','27','28'};

all_priority_id = [];
for jj = 1:length(runs);
    priority_id = [];
    filename = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/TASK/%s.mat',runs{jj})
    fileID = load(filename)
    priority_id = [priority_id; fileID.task.conditionAndQueriedTarget(:,1);];
    all_priority_id = [all_priority_id; priority_id]
    runnumtmp = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/DATA_PROC/mr_1_%s_proc.mat',runnum{jj});
    load(runnumtmp)
    
    ii_definetrial('XDAT',1,'XDAT',6)

for ii = 1:36 %:36 %[1 24];
tmp = [ii_cfg.trialvec XDAT]; %concatenate necessary vecs of equal length
test = tmp(:,1)==(ii) & tmp(:,2)==4;
XPri = X(test);
YPri = Y(test);
XPriend = XPri(end);
YPriend = YPri(end);
test2 = tmp(:,1)==(ii) & tmp(:,2)==5;
Xfinal = X(test2);
TarXnew = TarX(test2);
all_TarX = [all_TarX; TarXnew(1)];
TarYnew = TarY(test2);
all_TarY = [all_TarY; TarYnew(1)];


TarXnew = TarXnew(1);
TarYnew =TarYnew(1);

Xone = XPri(1);
Yone = YPri(1);

end 
end
%% consulted with tcs, try plotting keeping original circular mapping, final/ initial yoked to true loc
load('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/ii_stats')

all_TarX =[];
all_TarY =[];
all_Tar_rho =[];
%runs ={'run11','run12','run22','run23','run24','run25','run26','run27','run28'};
runs ={'run11'};

all_aligny = [];
all_alignx = [];
all_XPriend = []; 
all_YPriend = []; 

%runnum= {'11','12','22','23','24','25','26','27','28'};
runnum= {'11'}

all_priority_id = [];
for jj = 1:length(runs);
end
    priority_id = [];
    filename = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/TASK/%s.mat',runs{jj})
    fileID = load(filename)
    priority_id = [priority_id; fileID.task.conditionAndQueriedTarget(:,1);];
    all_priority_id = [all_priority_id; priority_id]
    runnumtmp = sprintf('/Volumes/hyper/experiments/Grace/TMS_Priority/subj03/l_ips2/DATA_PROC/mr_1_%s_proc.mat',runnum{jj});
    load(runnumtmp)
    
    ii_definetrial('XDAT',1,'XDAT',6)

 ii = 1%:36 %:36 %[1 24];

tmp = [ii_cfg.trialvec XDAT]; %concatenate necessary vecs of equal length
test = tmp(:,1)==(ii) & tmp(:,2)==4;
XPri = X(test);
Xpriend = XPri(end);

YPri = Y(test);
Ypriend = YPri(end);
test2 = tmp(:,1)==(ii) & tmp(:,2)==5;
Xfinal = X(test2);
TarXnew = TarX(test2);
all_TarX = [all_TarX; TarXnew(1)];
TarYnew = TarY(test2);
all_TarY = [all_TarY; TarYnew(1)];



TarXnew = TarXnew(1);
TarYnew =TarYnew(1);

priority_id_high = find(all_priority_id == 31);
high  = Xpriend(priority_id_high);
x_high = Xpriend(high);



      
priority_id_low = find(all_priority_id == 32);
low  = XPriend(priority_id_low);


figure(2);
plot(x_left_high,y_left_high,'b*','markersize',3)%,'markersize',4) %bo','markersize',4) hold on;
hold on;
plot(x_right_high,y_right_high,'go','markersize',3)%,'markersize',4) %bo','markersize',4)
plot(x_left_low,y_left_low,'r*','markersize',3)%,'markersize',4) %bo','markersize',4)
plot(x_right_low,y_right_low,'k*','markersize',3)%
hold on;
plot(avgx, avgy, 'k+','markersize',10)
