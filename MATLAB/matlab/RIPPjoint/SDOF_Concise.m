%   RIPP Joint Model, Segalman's Distribution, requires IMEX_a
clear
% Generic initialization of the joint properties
RIPP.Kp = 2e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = 4e3 ;
RIPP.KT = 1.5e7 ;
RIPP.chi = -0.5 ;
RIPP.beta = 0.005 ;
RIPP.theta = 1 ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;

% Hysteresis initialization
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 1 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   Hysteresis Loop Call
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nHysteresis Loop Calculations.\n')

x = 2.25e-3*sin(linspace(0,4*pi,72001));

tic
f2 = zeros(size(x)) ;
[f2(1), RIPP] = RIPPjoint(x(1),x(1)-0,RIPP) ;
for cntr = 2:length(x)
  [f2(cntr), RIPP] = RIPPjoint(x(cntr),x(cntr)-x(cntr-1),RIPP) ;
end
fprintf(['RIPP Segalman Distribution Time (red): ' t2str(toc) '.\n'])
figure(1)
plot(x,f2,'r')
xlabel('Displacement')
ylabel('Force')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   Dynamic System Integration Call
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nDynamic Simulations.\n')

% Define a sdof oscillator with no connections except the joint
% AMatrix is the state space representation of the system of equations
M = 1;
AMatrix = [0 0; 1 0] ;
Fappl = @(t) 5e3*sin(100*2*pi*t) ;

% Integration parameters
tstart = 0 ;
iparams.dt = 1e-6 ;
tfinal = 0.05 ;
iparams.tol = 1e-12 ;
iparams.dtmin = 1e-12 ;
iparams.dtmax = 1e-4 ;
iparams.output_res = 1e-5 ;
iparams.output_flag = 0 ;
iparams.cutoff = 600 ;
iparams.method = 5 ;

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 1 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Segalman Distribution Time (red): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),'r')
hold on
xlabel('Time')
ylabel('Displacement')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'r')
hold on
xlabel('Displacement')
ylabel('Force')

