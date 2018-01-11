addpath(genpath('/Volumes/hyper/experiments/Grace/iEye_irresponsible-master'))

%% make euclidean distance plot
load('/Volumes/hyper/experiments/Grace/tms_sessions/subj01/sham/DATA_PROC/1_5_proc.mat')
open iEye;
ii_definetrial

%%
runs ={'run05'};
fileweneed = [];
for jj = 1:length(runs);
filename = sprintf('/Volumes/hyper/experiments/Grace/tms_sessions/subj01/sham/TASK/%s.mat',runs{jj})
fileID = load(filename)
if  jj == 1;
    newrow = 1;
    endrow = 36;
else
newrow = ((jj-1).*36)+1
endrow = newrow + 35
end; 
fileweneed(newrow:endrow,1) = fileID.task.conditionAndQueriedTarget(:,1);
end 



for ii = 8
a = [ii_cfg.trialvec XDAT]; %concatenate necessary vecs of equal length
test = a(:,1)==(ii) & a(:,2)== 4;
Xnew = X(test);
Ynew = Y(test);
test2 = a(:,1)==(ii) & a(:,2)==5;
TarXnew = TarX(test2);
TarYnew = TarY(test2);
Xfinal = X(test2);
Yfinal = Y(test2);
TarXnew = unique(TarXnew(TarXnew~=0));
TarYnew = unique(TarYnew(TarYnew~=0));
 zpri = sqrt(Xnew.^2 +Ynew.^2);
 
if fileweneed(ii) == 31;  
 
% zfi = sqrt(Xfinal.^2 +Yfinal.^2);
%zt = sqrt(Xt.^2 + Yt.^2);
figure(1); 

plot(zpri,'b', 'linewidth', 2)
%plot(zfi,'r','linewidth',2)
else
    plot(zpri,'r', 'linewidth', 2)
end 
hold on; 
xlabel('Time(ms)')
ylabel('DVA')
ylim([0 11])
end 

% if fileweneed(ii) == 31; 
%     figure(1);
%     plot(0,0,'k+','markersize',10)
%     rxntime = ii_stats(1).srt(ii,1)/1000 %time in seconds
%     pause(rxntime)
%     plot(Xnew,Ynew,'b','linewidth', 1)
%     hold on;
% pause(0.1)
% plot(Xfinal, Yfinal,'b--')
% pause(0.1)
% plot(TarXnew,TarYnew,'ko','markersize',5,'markerfacecolor','b')
% else
%  
%     plot(Xnew, Ynew, 'r','linewidth', 1)
%     pause(0.1)
% plot(Xfinal, Yfinal,'r--')
% pause(0.1)
% plot(TarXnew,TarYnew,'ko','markersize',5,'markerfacecolor','r')
% end
% grid on
% xlim([-12 12])
% ylim([-12 12])
% xlabel('Horizontal DVA')
% ylabel('Vertical DVA')
% end 
% 
% z = sqrt(X.^2 + Y.^2)
