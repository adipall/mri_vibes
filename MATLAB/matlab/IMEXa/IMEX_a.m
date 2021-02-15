function [eta, t, F, UD] = IMEX_a(func, Amat, t0, t1, params)
% function [eta, t, F, UD] = IMEX_a(func, Amat, t0, t1, params)
% IMEX integration scheme using an adaptive time stepping algorithm.
%   calls RK1, RK2, RK3, RK4, RK5
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
%         method = 1 corresponds to the 1st order explicit Runge-Kutta
%            method termed the Euler method.  This uses a second order
%            method (the Midpoint rule) to adaptively update the time step.
%         method = 2 corresponds to the 2nd order explicit Runge-Kutta
%            method termed the Midpoint rule.  This uses a first order
%            method (Euler's method) to adaptively update the time step.
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
%                    Amat, tstart, tfinal, params) ;

%#ok<*CTCH>
if nargin < 4
  error('Too few input arguments.')
end
if nargin == 5
  if isfield(params,'dt')
    dt = params.dt ;
  else
    dt = (t1-t0)*1e-6 ;
  end
  if isfield(params,'dtmin')
    dtmin = params.dtmin ;
  else
    dtmin = 1e-12 ;
    if nargin < 5
      dt = max(dt,dtmin) ;
    else
      dtmin = min(dt,dtmin) ;
    end
  end
  if isfield(params,'dtmax')
    dtmax = params.dtmax ;
  else
    dtmax = 0.125*2*pi*max(abs(eigs(Amat))) ;
  end
  if isfield(params,'tol')
    tol = params.tol ;
  else
    tol = 1e-3 ;
  end
  if isfield(params,'IC')
    IC = params.IC ;
  else
    IC = zeros(length(Amat),1) ;
  end
  if isfield(params,'output_res')
    output_res = params.output_res ;
  else
    output_res = (t1-t0)/1000 ;
  end
  if isfield(params,'output_flag')
    output_flag = params.output_flag ;
  else
    output_flag = 1 ;
  end
  if isfield(params,'UO_func')
    UO_func = @(T,dt,Y,UD,flag) params.UO_func(T,dt,Y,UD,flag) ;
    out_flag = 1 ;
    clear UD
    UD.init = 1 ;
  else
    UD = 0 ;
    out_flag = 0 ;
  end
  if isfield(params,'ES_func')
    ES_func = @(T,Y,UD) params.ES_func(T,Y,UD) ;
    ES_flag = 1 ;
  else
    ES_flag = 0 ;
  end
  if isfield(params,'cutoff')
    cutoff = params.cutoff ;
  else
    cutoff = 3e6 ;
  end
  if isfield(params,'method')
    if params.method == 5
    % Uses RK54 Dormand-Prince method to compute next step, with adaptive
    % time step RK5(func, dt, dtmin, dtmax, tol, y, t)
      RKM = @(T,Y,UD,dt) RK5(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
    elseif params.method == 4
    % Uses RK45 Dormand-Prince method to compute next step, with adaptive
    % time step RK4(func, dt, dtmin, dtmax, tol, y, t)
      RKM = @(T,Y,UD,dt) RK4(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
    elseif params.method == 3
    % Uses RK32 Bogacki-Shampine method to compute next step, with adaptive
    % time step RK3(func, dt, dtmin, dtmax, tol, y, t)
      RKM = @(T,Y,UD,dt) RK3(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
    elseif params.method == 2
    % Uses RK21 Midpoint-Euler method to compute next step, with adaptive
    % time step RK1(func, dt, dtmin, dtmax, tol, y, t)
      RKM = @(T,Y,UD,dt) RK2(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
    elseif params.method == 1
    % Uses RK12 Euler-Midpoint method to compute next step, with adaptive
    % time step RK1(func, dt, dtmin, dtmax, tol, y, t)
      RKM = @(T,Y,UD,dt) RK1(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
    end
  else
    RKM = @(T,Y,UD,dt) RK5(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
  end
else
  dt = (t1-t0)*1e-6 ;
  dtmin = 1e-12 ;
  if nargin < 5
    dt = max(dt,dtmin) ;
  else
    dtmin = min(dt,dtmin) ;
  end
  dtmax = 1e-3 ;
  tol = 1e-3 ;
  IC = zeros(length(Amat),1) ;
  output_res = (t1-t0)/1000 ;
  output_flag = 1 ;
  UD = 0 ;
  out_flag = 0 ;
  ES_flag = 0 ;
  cutoff = 3e6 ;
  RKM = @(T,Y,UD,dt) RK5(@(t,y) func(t,y,UD), dt, dtmin, dtmax, tol, Y, T);
end

T = t0 ;
t_out = output_res+T ;
eta = zeros(ceil((t1-t0)/output_res),length(Amat)) ;
eta(1,:) = IC' ;
F = zeros(ceil((t1-t0)/output_res),length(Amat)) ;
t = zeros(ceil((t1-t0)/output_res),1) ;
S_flag = 0 ;
c_flag = 0 ;
Y = IC ;
cntr = 1 ;
eta(cntr,:) = Y' ;
t(cntr) = T ;
if isfield(UD,'init')
  UD = UO_func(T,dt,Y,UD,1) ;
  UD = rmfield(UD,'init') ;
end
F(cntr,:) = (func(T,Y,UD))' ;

ct1 = cputime ;

while T < t1 && S_flag == 0 ;
  % IMEX splitting method:

  % Forward explicit step
  [Y, T, dt1, err] = RKM(T,Y,UD,dt) ;

  % Backward Euler Method after the explicit step:
  Y = (eye(size(Amat))-dt*Amat) \ Y ;

  % Adjust dt for next step
  dt = dt1 ;

% Record output variables
  if T >= t_out
    t_out = t_out+output_res ;
    cntr = cntr+1 ;
    eta(cntr,:) = Y' ;
    t(cntr) = T ;
    F(cntr,:) = (func(T,Y,UD))' ;
    c_flag = 1 ;
    if out_flag == 1
      UD = UO_func(T,dt,Y,UD,1) ;
    end
  elseif out_flag == 1
    UD = UO_func(T,dt,Y,UD,0) ;
  end

% Output progress
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
    S_flag = ES_func(T,Y,UD) ;
  end

  if err > tol
    S_flag = 1 ;
  end

  if cputime > ct1+cutoff
    S_flag = 1 ;
  end
end

eta = eta(1:cntr,:) ;
t = t(1:cntr) ;
F = F(1:cntr,:) ;

