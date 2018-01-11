function [parspec] = myColors()
%for ease of getting the color scheme you want from bar plots 

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
% c11 = par(39,:,:)
% c12 = par(40,:,:)
c13 = par(55,:,:);
c14 = par(56,:,:);
c15 = brighten(c13,beta);
c16 = brighten(c14,beta);
% c15 = par(63,:,:)
% c16 = par(64,:,:)

parspec = [c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16];

end

