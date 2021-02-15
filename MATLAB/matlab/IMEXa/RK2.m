function [Y, T, dt, err] = RK2(func, dt, dtmin, dtmax, tol, y, t)
% Second order Runge Kutta explicit step, for use in an IMEX algorithm.
% Uses an Euler first order and midpoint second order scheme.
% Adaptive time stepping provided by the error between the two methods.

i = 1 ;
s = 1 ;

% Start Iterations
while i<2
  % Modify Time Step
  dt = dt*s ;

  % 1st order method
  F1 = func(t+dt,y) ;
  y1 = y+dt*F1 ;
  % 2nd order function call
  F2 = func(t+dt/2,y+1/2*dt*F1) ;
  z1 = y+dt*F2 ;

  err = abs(z1-y1) ;

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

Y=z1;

clear i j k1 k2 k3 k4 k5 k6 k7 s

