%   Original Segalman Iwan Joint Models, requires IMEX_a
% number of joint degrees of freedom: Njoints
joint_states = zeros(30,1) ; % the states of the jenkins elements
joint_params = zeros(1,4) ; %#ok<PREALL>

% generic initialization of the joint parameters
% JointsData contains the matrix YL, which is size 4x1994
% load JointsData
% joint_Number = round(rand*1993+1);
% joint_params(1,:) = YL(:,joint_Number) ;
% % params = [ chi, beta, K_T, F_S]
joint_params = [-0.5 0.005 1.5e7 4e3] ;

params.joint = joint_params ;
params.y0 = joint_states ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Segalman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generic initialization of the joint properties
RIPP.Kp = 0;%1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.theta = 1;

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
f = zeros(size(x)) ;
for cntr = 1:length(x)
  [f(cntr), params] = IMEX_Iwan(x(cntr),params);
end
fprintf(['Segalman''s Iwan Model Time (blue): ' t2str(toc) '.\n'])

figure(1)
plot(x,f,'b')
hold on
xlabel('Displacement')
ylabel('Force')

tic
f2 = zeros(size(x)) ;
[f2(1), RIPP] = RIPPjoint(x(1),x(1)-0,RIPP) ;
for cntr = 2:length(x)
  [f2(cntr), RIPP] = RIPPjoint(x(cntr),x(cntr)-x(cntr-1),RIPP) ;
end
fprintf(['RIPP Segalman Distribution Time (red): ' t2str(toc) '.\n'])
figure(1)
plot(x,f2,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Mignolet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.Kp = 1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
RIPP.model = 1 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

tic
f3 = zeros(size(x)) ;
[f3(1), RIPP] = RIPPjoint(x(1),x(1)-0, RIPP) ;
for cntr = 2:length(x)
  [f3(cntr), RIPP] = RIPPjoint(x(cntr),x(cntr)-x(cntr-1),RIPP) ;
end
fprintf(['RIPP Mignolet Distribution Time (green): ' t2str(toc) '.\n'])
figure(1)
plot(x,f3,'g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Uniform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.Kp = 1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 2 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

tic
f4 = zeros(size(x)) ;
[f4(1), RIPP] = RIPPjoint(x(1),x(1)-0,RIPP) ;
for cntr = 2:length(x)
  [f4(cntr), RIPP] = RIPPjoint(x(cntr),x(cntr)-x(cntr-1),RIPP) ;
end
fprintf(['RIPP Uniform Distribution Time (magenta): ' t2str(toc) '.\n'])
figure(1)
plot(x,f4,'m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Iwan-Stribeck, v scale = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.Kp = 1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
Fv = 1e-3 ;
vs = 1e-4 ;
ds = 2 ;
RIPP.Fv = Fv ;
RIPP.ds = ds ;
RIPP.vs = vs ;
RIPP.model = 3 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

tic
f5 = zeros(size(x)) ;
[f5(1), RIPP] = RIPPjoint(x(1),(x(1)-0)*1, RIPP) ;
for cntr = 2:length(x)
  [f5(cntr), RIPP] = RIPPjoint(x(cntr),(x(cntr)-x(cntr-1))*1,RIPP) ;
end
fprintf(['RIPP Iwan-Stribeck (v scale = 1) Time (black solid): ' t2str(toc) '.\n'])
figure(1)
plot(x,f5,'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Iwan-Stribeck, v = 1e2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.Kp = 1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
RIPP.Fv = Fv ;
RIPP.ds = ds ;
RIPP.vs = vs ;
RIPP.model = 3 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

tic
f6 = zeros(size(x)) ;
[f6(1), RIPP] = RIPPjoint(x(1),(x(1)-0)*1e2, RIPP) ;
for cntr = 2:length(x)
  [f6(cntr), RIPP] = RIPPjoint(x(cntr),(x(cntr)-x(cntr-1))*1e2,RIPP) ;
end
fprintf(['RIPP Iwan-Stribeck (v scale = 1e2) Time (black dashed): ' t2str(toc) '.\n'])
figure(1)
plot(x,f6,'--k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Iwan-Stribeck, v = 1e4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.Kp = 1e7 ;
RIPP.dp = 2e-3 ;
RIPP.FS = joint_params(4) ;
RIPP.KT = joint_params(3) ;
RIPP.chi = joint_params(1) ;
RIPP.beta = joint_params(2) ;
RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
RIPP.Fv = Fv ;
RIPP.ds = ds ;
RIPP.vs = vs ;
RIPP.model = 3 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

tic
f7 = zeros(size(x)) ;
[f7(1), RIPP] = RIPPjoint(x(1),(x(1)-0)*1e4, RIPP) ;
for cntr = 2:length(x)
  [f7(cntr), RIPP] = RIPPjoint(x(cntr),(x(cntr)-x(cntr-1))*1e4,RIPP) ;
end
fprintf(['RIPP Iwan-Stribeck (v scale = 1e4) Time (black dotted): ' t2str(toc) '.\n'])
figure(1)
plot(x,f7,':k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   Dynamic System Integration Call
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nDynamic Simulations.\n')

M = 1;

AMatrix = [0 0; 1 0] ;
Fappl = @(t) 5e3*sin(15*2*pi*t) ;
Force_Func = @(t,y,UD) [(Fappl(t)-IMEX_Iwan(y(2),UD))/M; 0] ;

tstart = 0 ;
iparams.dt = 1e-6 ;
tfinal = 0.3 ;
iparams.tol = 1e-12 ;
iparams.dtmin = 1e-12 ;
iparams.dtmax = 1e-4 ;
iparams.output_res = 1e-5 ;
iparams.output_flag = 0 ;
iparams.cutoff = 600 ;
iparams.method = 5 ;

params.y0 = joint_states ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_Iwan(y(2),UD,params) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['Segalman''s Iwan Model Time (blue): ' t2str(toc) '.\n'])

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   Original Segalman Iwan Joint Models
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

joint_function_handle = @iwan_function_preferred ;
% number of joint degrees of freedom: Njoints
joint_states = zeros(500,1) ;

params.joint = joint_params ;
params.y0 = joint_states ;

Force_Func = @(t,y,UD) [(Fappl(t)-IMEX_Iwan(y(2),UD))/M; 0] ;

iparams.UO_func = @(t,dt,y,UD,flag) UO_Iwan(y(2),UD,params) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['Segalman''s Iwan Model Time (500 elements, cyan): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),'c')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'c')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'c')

figure(4)
plot(Y(:,2),Y(:,1),'c')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'c')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Segalman
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 1 ;
RIPP.theta = 1 ;
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

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'r')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'r')

figure(4)
plot(Y(:,2),Y(:,1),'r')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Mignolet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 1 ;
RIPP.theta = theta ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;


Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Mignolet Distribution Time (green): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),'g')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'g')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'g')

figure(4)
plot(Y(:,2),Y(:,1),'g')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'g')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Uniform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
RIPP.model = 2 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Uniform Distribution Time (magenta): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),'m')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'m')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'m')

figure(4)
plot(Y(:,2),Y(:,1),'m')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Iwan-Stribeck
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
RIPP.model = 3 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;

Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Iwan-Stribeck (black): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),'k')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),'k')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),'k')

figure(4)
plot(Y(:,2),Y(:,1),'k')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   RIPP Joint Model, Iwan-Stribeck
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RIPP.F0 = 0 ;
RIPP.d0 = 0 ;
theta = 0.75 ; % Mignolet's fifth parameter
RIPP.theta = theta ;
RIPP.model = 3 ;
RIPP.direction = 0 ;
RIPP.Flast = 0 ;
RIPP.dlast = 0 ;
Fv = 100 ;
vs = 1e-4 ;
ds = 100 ;
RIPP.Fv = Fv ;
RIPP.ds = ds ;
RIPP.vs = vs ;

Force_Func = @(t,y,UD) [(Fappl(t)-RIPPjoint(y(2),y(1),UD))/M; 0] ;
iparams.UO_func = @(t,dt,y,UD,flag) UO_RIPP(y(2),y(1),UD,RIPP) ;

tic
[Y, T, F, ~] = IMEX_a(@(t,y,UD) Force_Func(t,y,UD), ...
                       AMatrix, tstart, tfinal, iparams) ;
fprintf(['RIPP Iwan-Stribeck (black): ' t2str(toc) '.\n'])

figure(2)
plot(T,Y(:,2),':k')

figure(3)
plot(Y(:,2),-(F(:,1)-Fappl(T)),':k')

figure(6)
plot(Y(round(end/2):end,2),-(F(round(end/2):end,1)-Fappl(T(round(end/2):end))),':k')

figure(4)
plot(Y(:,2),Y(:,1),':k')

figure(5)
plot(Y(round(end/2):end,2),Y(round(end/2):end,1),':k')

fprintf('\n')
