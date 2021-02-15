function [eta, eta2, t, F, UD] = IMEX_2a(func, M, C, K, t0, t1, params)
import imex_a.*;   % TODO figure out what this means
% IMEX integration scheme using an adaptive time stepping algorithm.
%    calls RKNG4 and RKNG5
% This scheme uses a fifth order Runge-Kutta method for the explicit step,
% a backward Euler method for the implitic step, and an RK45 Fehlberg
% method to determine the time step size.
%
% MR Brake, IMEXa: An adaptive 5th order implicit-explicit integration
% scheme, sand2013-4299, 2013
%
% This function takes 5 input arguments:
% func : the force function for the simulation
% Amat : the system matrix for the simulation
% t0 : the starting time
% t1 : the final time
% params : a structure containing the following parameters:
%   dt : the initial time step
%   dtmin : the minimum time step size
%   dtmax : the maximum time step size
%   tol : the tolerance for the error in calculating the adaptive time step
%   IC : the initial conditions for the temporal coefficients for the
%        degreesof freedom
%   output_res : the temporal resolution of the output
%   output_flag : determines whether or not to display output steps
%   UO_func : the user defined output function that specifies quantities to
%             be saved into the UD output structure.
%   ES_func : an early stop criterion for the simulation
%   cutoff : in case of poor convergence, this is a specified cpu wall time
%            to terminate the simulation at.
%   method : explicit/adaptive time stepping method
%         method = 3 corresponds to the 3rd order explicit Runge-Kutta
%            method termed the Bogacki-Shampine method.  This uses a second
%            order approximation to adaptively update the time step.
%         method = 5 (the default) corresponds to the 5th order explicit
%            Runge-Kutta method termed the Dormand-Prince method.  The time
%            step is adaptively chosen using a fourth order estimation.
%
% The first four inputs are strictly required.
% ES_func and UO_func, if not supplied, are disabled in the simulation.
% All other parameters are given estimated values.  These values will not
% be optimal for a generic simulation, and the user is strongly encouraged
% to determine appropriate values for each input.
% If UO_func is not specified, UD is returned as 0.
%
% Output for this function are divided into four structures:
% eta : the time history of the temporal coefficients for the degrees of
%       freedom
% t : the time history record
% F : the force, calculated by the user defined force function
% UD : a structure calculated by the user defined output function.
%
% An example of calling this function is:
%
% [eta, t, F, UD] = IMEX_a(@(t,y,UD) my_force_func(t,y,UD,params), ...
%                    Amat, tstart, tfinal, dt, dtmin, dtmax, ats_tol, ...
%                    IC, t_out, output_flag, ...
%                    @(t,y,UD,flag) my_Output_func(t,dt,y,UD,flag,params), ...
%                    @(t,y,UD) my_ES_func(t,y,UD,params)) ;

%#ok<*CTCH>
if nargin < 6
  error('Too few input arguments.')
end

% Set up system matrices
G = M\K;
H = M\C;

Force = @(t,y,dy,UD) func(t, y, dy, UD);

% Defaults
dt = (t1-t0)*1e-6 ;
dtmin = 1e-12 ;
dt = max(dt, dtmin) ;
dtmax = 1e-3 ;
tol = 1e-3 ;
IC = zeros(length(M),1) ;
IC2 = zeros(length(M),1) ;
output_res = (t1-t0)/1000 ;
output_flag = 1 ;
UD = 0 ;
out_flag = 0 ;
ES_flag = 0 ;
cutoff = 3e6 ;
RKx = @RKNG5 ;

if nargin == 7
  if isfield(params,'dt'), dt = params.dt; end
  if isfield(params,'dtmin')
    dtmin = params.dtmin;
  else
    dtmin = min(dt,dtmin) ;
  end
  if isfield(params,'dtmax'), dtmax = params.dtmax; end
  if isfield(params,'tol'), tol = params.tol; end
  if isfield(params,'IC'), IC = params.IC; end
  if isfield(params,'IC2'), IC2 = params.IC2; end
  if isfield(params,'output_res'), output_res = params.output_res; end
  if isfield(params,'output_flag'), output_flag = params.output_flag; end
  if isfield(params,'UO_func')
    UO_func = @(T,dt,Y,DY,UD,flag) params.UO_func(T,dt,Y,DY,UD,flag) ;
    out_flag = 1 ;
    clear UD
    UD.init = 1 ;
  end
  if isfield(params,'ES_func')
    ES_func = @(T,Y,DY,UD) params.ES_func(T,Y,DY,UD) ;
    ES_flag = 1 ;
  end
  if isfield(params,'cutoff'), cutoff = params.cutoff; end
  if isfield(params,'method')
    switch params.method
      case 5  % Uses RKNG45 method
        RKx = @RKNG5;
      case 4  % Uses RKNG34 method
        RKx = @RKNG4;
      otherwise
        RKx = @RKNG5;
    end
  end
end

RKM = @(T,Y,DY,UD,dt) RKx(@(t,y,dy) func(t,y,dy,UD), M, G, H, dt, dtmin, dtmax, tol, Y, DY, T);

T = t0 ;
t_out = output_res+T ;

eta = zeros(ceil((t1-t0)/output_res),length(M)) ;
eta2 = zeros(ceil((t1-t0)/output_res),length(M)) ;
F = zeros(ceil((t1-t0)/output_res),length(M)) ;
t = zeros(ceil((t1-t0)/output_res),1) ;

S_flag = 0 ;
c_flag = 0 ;
Y = IC ;
DY = IC2 ;
cntr = 1 ;

eta(cntr,:) = Y' ;
eta2(cntr,:) = DY' ;
t(cntr) = T ;
if isfield(UD,'init')
  UD = UO_func(T,dt,Y,DY,UD,1) ;
  UD = rmfield(UD,'init') ;
end
F(cntr,:) = (func(T,Y,DY,UD))' ;

ct1 = cputime ;

while T < t1 && S_flag == 0 ;
  % IMEX splitting method:
  oldY = Y;

  % Forward explicit step
  [allY, T, dt, dt_new, err] = RKM(T,Y,DY,UD,dt) ;
  Y = allY.Y ;
  DY = allY.DY ;

  %Implicit step; this is the Trapezoidal rule...
  T = T + dt;  % This accounts for the RKM returning Y(n+1/2)
  dt = 2 * dt;
  tstep = 1/(.5*dt);
  Y = (tstep^2*M + tstep*C + K)\( Force(T,Y,DY,UD) - tstep^2*M*(oldY - 2*Y) + tstep*C*Y);
  acc = M\(Force(T, Y, allY.DY,UD)) - H*allY.DY - G*Y;
  DY = allY.DY + acc*(dt/2);

  dt = dt_new;


% Record output variables
  if T >= t_out
    t_out = t_out+output_res ;
    cntr = cntr+1 ;
    eta(cntr,:) = Y' ;
    eta2(cntr,:) = DY' ;
    t(cntr) = T ;
    F(cntr,:) = (func(T,Y,DY,UD))' ;
    c_flag = 1 ;
    if out_flag == 1
      UD = UO_func(T,dt,Y,DY,UD,1) ;
    end
  elseif out_flag == 1
    UD = UO_func(T,dt,Y,DY,UD,0) ;
  end

% Output progres
  if output_flag == 1
    if cntr/2e1 == floor(cntr/2e1) && c_flag == 1
      fprintf('.')
      if cntr/1e3 == floor(cntr/1e3)
        fprintf('\n')
      end
      c_flag = 0 ;
    end
  end

% Check for early termination conditions
  if ES_flag == 1
    S_flag = ES_func(T,Y,DY,UD) ;
  end

  if norm(err) > tol && dt <= dtmin
    S_flag = 1 ;
  end

  if cputime > ct1+cutoff
    S_flag = 1 ;
  end
end

eta = eta(1:cntr,:) ;
eta2 = eta2(1:cntr,:) ;
t = t(1:cntr) ;
F = F(1:cntr,:) ;

end
