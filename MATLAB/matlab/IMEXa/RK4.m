function [Y, T, dt, err] = RK4(func, dt, dtmin, dtmax, tol, y, t)
% Fifth order Runge Kutta explicit step with an adaptive time step
% (Dormand-Prince method), for use in an IMEX algorithm.

i = 1;
s = 1 ;

% Start Iterations
while i~=2
  % Modify Time Step
  dt = dt*s ;

  % RK45 Dormand-Prince Method
  k1 = dt * func(t,y);
  k2 = dt * func(t+1/5*dt, y+(1/5)*k1);
  k3 = dt * func(t+3/10*dt, y+(3/40)*k1+(9/40)*k2);
  k4 = dt * func(t+4/5*dt, y+(44/45)*k1-(56/15)*k2+(32/9)*k3);
  k5 = dt * func(t+8/9*dt, y+(19372/6561)*k1-(25360/2187)*k2+(64448/6561)*k3-(212/729)*k4);
  k6 = dt * func(t+dt, y+(9017/3168)*k1-(355/33)*k2+(46732/5247)*k3+(49/176)*k4-(5103/18656)*k5);

  % Fifth order accurate solution:
  y1 = y+(35/384)*k1+(500/1113)*k3+(125/192)*k4-(2187/6784)*k5+(11/84)*k6;

  k7 = dt * func(t+dt,y1);
  % Fourth order accurate solution:
  Y = y + (5179/57600)*k1+(7571/16695)*k3+(393/640)*k4-(92097/339200)*k5+(187/2100)*k6+(1/40)*k7;

  % Error = Fifth order solution - Fourth order solution:
  err = abs((35/384-5179/57600)*k1+(500/1113-7571/16695)*k3+(125/192-393/640)*k4-(2187/6784-92097/339200)*k5+(11/84-187/2100)*k6+(0-1/40)*k7);

  % Time Step criteria, differs from Fehlberg from trial & error
  s=((tol)/(norm(err)))^.25;

  % Create a window for s
  if s<(1/2)^(1/4) % Error is twice tolerance
    j = 1 ; % Marker for redo-ing a step
    if dt<=dtmin
      warning(['Minimum time step too large, tolerance not met at time ' num2str(1000*t) ' ms.']) %#ok<WNTAG>
      j = 0 ;
      err = tol*2 ;
    end
  elseif s>sqrt(10)
    if sqrt(10)*dt >= dtmax
      s = dtmax/dt ;
      j = 0 ;
    else
      s = sqrt(10);
      j = 0;
    end
  elseif s<.9
    if dt <= dtmin
      j = 0;
    else
      j = 1;
    end
  else
    j=0;
  end

  if j==1
    if s < 0.85
      i = 1 ; % Repeats the last step with a different timestep
    else
      i = i+1 ;
    end
    % Place a window on the time step
    if dt < dtmin
      dt = dtmin ;
    elseif dt > dtmax
      dt = dtmax ;
    end
  else
    i = 2 ;
  end
end
T=t+dt;
dt = dt*s ;

clear i j k1 k2 k3 k4 k5 k6 k7 s

