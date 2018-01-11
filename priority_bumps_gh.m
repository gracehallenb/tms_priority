
[x,y]=meshgrid(-3:0.1:3, -3:0.1:3); 
%% four item priority
figure(1)
subplot(1,2,1)
surf(x, y, 0 + ((1*(x-(cosd(90))).^2+1.25)./(10*(x-(cosd(90))).^2+1.25)).*((1*(y-(sind(0)*2)).^2+1.25)./(10*(y-(sind(0)*2)).^2+1.25))...
     + ((-1*(x-(cosd(45)*2)).^2+2)./(10*(x-(cosd(45)*2)).^2+1)).*((-1*(y-(sind(45)*2)).^2+2)./(10*(y-(sind(45)*2)).^2+1))...  %2
 + ((-1*(x-(cosd(45)*2)).^2+1.3)./(10*(x-(cosd(45)*2)).^2+1)).*((-1*(y-(-sind(45)*2)).^2+1.3)./(10*(y-(-sind(45)*2)).^2+1))... %5
 + ((-1*(x-(-cosd(45)*2)).^2+1.75)./(10*(x-(-cosd(45)*2)).^2+1)).*((-1*(y-(-sind(45)*2)).^2+1.75)./(10*(y-(-sind(45)*2)).^2+1))... %8
 + ((-1*(x-(-cosd(45)*2)).^2+1)./(10*(x-(-cosd(45)*2)).^2+1)).*((-1*(y-(sind(45)*2)).^2+1)./(10*(y-(sind(45)*2)).^2+1))) %11

view(0,65)
title('Four Item Priority ')
shading interp
axis([-3 3 -3 3 -1 5]) %fixation

% two item priority

figure (1)
subplot(1,2,2)
surf(x, y, 0 + ((1*(x-(cosd(90))).^2+1.25)./(10*(x-(cosd(90))).^2+1.25)).*((1*(y-(sind(0)*2)).^2+1.25)./(10*(y-(sind(0)*2)).^2+1.25))...
     + ((-1*(x-(cosd(45)*2)).^2+2.25)./(10*(x-(cosd(45)*2)).^2+1)).*((-1*(y-(sind(45)*2)).^2+2.25)./(10*(y-(sind(45)*2)).^2+1))...  %2
+ ((-1*(x-(-cosd(45)*2)).^2+1.9)./(10*(x-(-cosd(45)*2)).^2+1)).*((-1*(y-(-sind(45)*2)).^2+1.9)./(10*(y-(-sind(45)*2)).^2+1)))%8

view(0,65)
title('Two Item Priority ')
axis([-3 3 -3 3 -1 5]) %fixation

%% with tms 


figure (2)

% this plot reflects two item representation before tms
% subplot(1,3,1)
% h = surf(x, y, 0 + ((1*(x-(cosd(90))).^2+0.5)./(10*(x-(cosd(90))).^2+0.5)).*((1*(y-(sind(0)*2)).^2+0.5)./(10*(y-(sind(0)*2)).^2+0.5))...
%      + ((-1*(x-(cosd(45)*2)).^2+2.25)./(10*(x-(cosd(45)*2)).^2+1)).*((-1*(y-(sind(45)*2)).^2+2.25)./(10*(y-(sind(45)*2)).^2+1))...  %2
% + ((-1*(x-(-cosd(45)*2)).^2+1.8)./(10*(x-(-cosd(45)*2)).^2+1)).*((-1*(y-(-sind(45)*2)).^2+1.8)./(10*(y-(-sind(45)*2)).^2+1)))%8
% %set(h,'linestyle','none')
% %shading interp
% view(0,65)
% title('Two Item Priority, no TMS')
% axis([-3 3 -3 3 -1 5]) %fixation

% this plot reflects two item representation aftertms -- need to change variance here
subplot(1,2,1)
g = surf(x, y, 0 + ((1*(x-(cosd(90))).^2+0.5)./(20*(x-(cosd(90))).^2+1)).*((1*(y-(sind(0)*2)).^2+0.5)./(10*(y-(sind(0)*2)).^2+ 1))...
     + ((-1*(x-(cosd(45)*2)).^2+2.025)./(10*(x-(cosd(45)*2)).^2+2.025)).*((-1*(y-(sind(45)*2)).^2+2.025)./(10*(y-(sind(45)*2)).^2+2.025))...  %2
+ ((-1*(x-(-cosd(45)*2)).^2+2.025)./(10*(x-(-cosd(45)*2)).^2+2.025)).*((-1*(y-(-sind(45)*2)).^2+2.025)./(10*(y-(-sind(45)*2)).^2+2.025)))%8
%set(g,'linestyle','none')
%shading interp
view(0,65)
title('Two Item Priority, sPCS TMS ')
axis([-3 3 -3 3 -1 5]) %fixation


subplot(1,2,2)
g = surf(x, y, 0 + ((1*(x-(cosd(90))).^2+0.5)./(20*(x-(cosd(90))).^2+1)).*((1*(y-(sind(0)*2)).^2+0.5)./(10*(y-(sind(0)*2)).^2+1))...
     + ((-1*(x-(cosd(45)*2)).^4+2.25)./(10*(x-(cosd(45)*2)).^4+2)).*((-1*(y-(sind(45)*2)).^4+ 2.25)./(10*(y-(sind(45)*2)).^4+2))...  %2
+ ((-1*(x-(-cosd(45)*2)).^4+1.8)./(10*(x-(-cosd(45)*2)).^4+1.8)).*((-1*(y-(-sind(45)*2)).^4+1.8)./(10*(y-(-sind(45)*2)).^4+1.8)))%8
%set(g,'linestyle','none')
%shading interp
view(0,65)
title('Two Item Priority, sPCS TMS ')
axis([-3 3 -3 3 -1 5]) %fixation

