function [Y, T, dt, dt2, err] = RKNG5(func, M, G, H, dt, dtmin, dtmax, tol, y, dy, t)
% RKNG 45.  See Fine and Haute (1985).

level = 7;

a = [[                    0,                0,                  0,                  0,                 0, 0, 0];
            [              32/1521,                0,                  0,                  0,                 0, 0, 0];
            [                4/169,            4/169,                  0,                  0,                 0, 0, 0];
            [             175/5184,                0,          1625/5184,                  0,                 0, 0, 0];
            [-342497279/5618900760, 6827067/46824173, 35048741/102161832, -2201514/234120865,                 0, 0, 0];
            [          -7097/52152,         767/2173,        14027/52152,            30/2173,                 0, 0, 0];
            [           4817/51600,                0,     388869/1216880,         3276/23575, -1142053/22015140, 0, 0]];

da = [[                 0,                   0,                  0,                   0,                  0,      0, 0];
             [              8/39,                   0,                  0,                   0,                  0,      0, 0];
             [              1/13,                3/13,                  0,                   0,                  0,      0, 0];
             [         7385/6912,          -9425/2304,         13325/3456,                   0,                  0,      0, 0];
             [223324757/91364240, -174255393/18272848, 382840094/46824173, -39627252/234120865,                  0,      0, 0];
             [      108475/36464,           -9633/848,     7624604/806183,          8100/49979,  -4568212/19446707,      0, 0];
             [        4817/51600,                   0,    1685099/3650640,         19656/23575, -53676491/88060560, 53/240, 0]];

b     = ([  4817/51600,    0,     388869/1216880,       3276/23575,     -1142053/22015140,          0,     0])';
db    = ([  4817/51600,    0,    1685099/3650640,      19656/23575,    -53676491/88060560,     53/240,     0])';
bHat  = ([8151/2633750,    0, -1377519/186334750,  586872/28879375,  -36011118/2247378875,          0,     0])';
dbHat = ([8151/2633750,    0, -5969249/559004250, 3521232/28879375, -846261273/4494757750, 4187/36750, -1/25])';
c     =  [           0, 8/39,               4/13,              5/6,                 43/47,          1,     1];


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
