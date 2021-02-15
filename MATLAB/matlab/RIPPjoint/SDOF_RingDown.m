%   RIPP Joint Model, Segalman's Distribution, requires IMEX_a
% Generic initialization of the joint properties
RIPP.Kp = 2e7 ;
RIPP.dp = 2e-3 ;
% RIPP.FS = 4e3 ;
% RIPP.KT = 1.5e7 ;
% RIPP.chi = -0.5 ;
% RIPP.beta = 0.005 ;
RIPP.FS = 50 ;
RIPP.KT = 1.5e6 ;
RIPP.chi = -0.05 ;
RIPP.beta = 0.01 ;
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
%%   Dynamic System Integration Call
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nDynamic Simulations.\n')

% Define a sdof oscillator with no connections except the joint
% AMatrix is the state space representation of the system of equations
M = 1;
AMatrix = [0 0; 1 0] ;
Fappl = @(t) 0;% 5e3*(1-(t >= 0.01)) ;
iparams.IC = [0.02; 0];

% Integration parameters
tstart = 0 ;
iparams.dt = 1e-6 ;
tfinal = 0.5 ;
iparams.tol = 1e-12 ;
iparams.dtmin = 1e-16 ;
iparams.dtmax = 4e-6 ;
iparams.output_res = 1e-6 ;
iparams.output_flag = 1 ;
iparams.cutoff = 6000 ;
iparams.method = 5 ;

Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Segalman Distribution Time: ' t2str(toc) '.\n'])


figure(2)
plot(T,Y(:,2),'b')
hold on
xlabel('Time')
ylabel('Displacement')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'b')
hold on
xlabel('Displacement')
ylabel('Force')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'b')
hold on
xlabel('Displacement')
ylabel('Force')

figure(4)
plot(Y(:,2),Y(:,1),'b')
hold on
xlabel('Displacement')
ylabel('Velocity')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'b')
hold on
xlabel('Displacement')
ylabel('Velocity')

Ti = linspace(0,tfinal,floor(tfinal/iparams.output_res)+1) ;
Yi = interp1(T,Y(:,2),Ti,'spline');
