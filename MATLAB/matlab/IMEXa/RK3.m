function [Y, T, dt, err] = RK3(func, dt, dtmin, dtmax, tol, y, t)
% Third order Runge Kutta explicit step with an adaptive time step
% (Bogacki-Shampine method), for use in an IMEX algorithm.

i = 1 ;
s = 1 ;

% Start Iterations
while i<2
  % Modify Time Step
  dt = dt*s ;

  % RK3 Bogacki-Shaampine Method
  k1 = func(t,y);
  k2 = func(t+1/2*dt, y+(1/2)*dt*k1);
  k3 = func(t+3/4*dt, y+(3/4)*dt*k2);
  y1 = y + 2/9*dt*k1+1/3*dt*k2+4/9*dt*k3 ;
  k4 = func(t+dt, y1);
  z1 = y + 7/24*dt*k1+1/4*dt*k2+1/3*dt*k3+1/8*dt*k4 ;

  err = abs(z1-y1);

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

Y=y1;

clear i j k1 k2 k3 k4 k5 k6 k7 s

