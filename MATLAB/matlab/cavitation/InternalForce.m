%% ------------------------------------------------------------------------
% This is going to be a script to test out the Internal Force calculation
% for cavitating elements
%
% MRR, 1/3/2017
%% ------------------------------------------------------------------------

% ---- Input parameters

% Coordinates
ex = [0,1,1,0,0,1,1,0];
ey = [0,0,1,1,0,0,1,1];
ez = [0,0,0,0,1,1,1,1];

rho = 1.23;cd
c = 2.34;
pvapor = 0.2;
penalty = 3.45;

vel_pot = ex;
dvel_pot = -ey;
ddvel_pot = ez;

ep = 2; % integration rule

%% ---- Here we go
[Mf,mass] = soli8ma(ex,ey,ez,ep,rho,c); %% actually wrong... shouldn't it be 1/c^2*tM instead of rho*tM?
if (mass == rho)
    disp('Brilliant');
else
    disp('Blame Manoj');
end
Kf = soli8Ka(ex,ey,ez,ep,rho,c);
Bf = soli8Ba(ex,ey,ez,ep,pvapor,dvel_pot,penalty);

Internal_force = (1/c^2)*Mf*ddvel_pot' + Bf*(dvel_pot-pvapor)' + Kf*vel_pot';

tv1=Mf*ddvel_pot';
tv2=Bf*(dvel_pot-pvapor)';
tv3=Kf*vel_pot';

tv01=1/(24*c^2)*[1;1;1;1;2;2;2;2];
tv02=-penalty*(1/24*([1;1;2;2;1;1;2;2])+pvapor/8);
tv03=1/4*[-1;1;1;-1;-1;1;1;-1];
