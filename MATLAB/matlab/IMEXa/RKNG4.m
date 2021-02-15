function [Y, T, dt, dt2, err] = RKNG4(func, M, G, H, dt, dtmin, dtmax, tol, y, dy, t)
% RKNG 34.  See Fine and Haute (1985).

level = 5;

a = [[    0,    0,       0,    0, 0];
            [ 2/81,    0,       0,    0, 0];
            [ 1/36,  1/36,      0,    0, 0];
            [9/128,     0, 27/128,    0, 0];
            [11/60, -3/20,   9/25, 8/75, 0]];

da = [[     0,        0,      0,     0, 0];
             [   2/9,        0,      0,     0, 0];
             [  1/12,      1/4,      0,     0, 0];
             [69/128, -243/128, 135/64,     0, 0];
             [-17/12,     27/4,  -27/5, 16/15, 0]];

b     = ([ 19/180,   0,   63/200,  16/225,   1/120])';
db    = ([    1/9,   0,     9/20,   16/45,    1/12])';
bHat  = ([25/1116,   0, -63/1240, 64/1395, -13/744])';
dbHat = ([  2/125,   0,  -27/625,  32/625,  -3/125])';
c     =  [      0, 2/9,      1/3,     3/4,       1];


redo = 1;
s = 1;

while redo
  dt = dt*s ;

  f = zeros(length(y), floor(level));
  for i=1:level
    sumAf = f * a(i,:)';
    sumDaF = f * da(i,:)';

    f(:,i) = M\(func(t+c(i)*dt, y+c(i)*dt*dy + dt^2*sumAf, dy + dt*sumDaF)) - H*(dy + dt*sumDaF) - G*(y+c(i)*dt*dy + dt^2*sumAf);
  end

  Y.Y = y + dt*dy + dt^2*(f*b);
  Y.DY = dy + dt*(f*db);

  err = abs(dt^2*f*bHat) + abs(dt*f*dbHat);


  % Time Step criteria, differs from Fehlberg from trial & error
  ne = norm(err);
  if ne - 0 < eps
      s = 10;
  else
      s=(tol/ne)^.25;
  end

% Create a window for s
  if s<(1/2)^(1/4) % Error is twice tolerance
    redo = 1 ; % Marker for redo-ing a step
    if dt<=dtmin
      warning(['Minimum time step too large, tolerance not met at time ' num2str(1000*t) ' ms.']) %#ok<WNTAG>
      redo = 0 ;
      err = tol*2 ;
    end
  elseif s>sqrt(10)
    redo = 0;
    if sqrt(10)*dt >= dtmax
      s = dtmax/dt ;
    else
      s = sqrt(10);
    end
  else
    redo=0;
  end

  if redo
    % Place a window on the time step
    if dt < dtmin
      dt = dtmin ;
    elseif dt > dtmax
      dt = dtmax ;
    end
  end
end
T=t+dt;
dt2 = dt*s ;


end
