%TESTCONSTANTFORCE test solutionConstantForce
clear
[k,c,m] = getRandomModel();
uforce = randn();
n=1000;
tfinal = 1;
[u,vel] = solutionConstantForce(m,c,k,uforce,n,tfinal);
time = linspace(0,tfinal,n);
figure(1); plot(time,u);title('disp');
dt = time(2);
velocity = (u(2:end) - u(1:end-1))/dt;
t2 = .5*(time(2:end) + time(1:end-1));
figure(2);
plot(t2,velocity,'-.', time ,vel,'.');
title('velocity');
v2 = .5*(velocity(2:end) + velocity(1:end-1));
acceleration =(velocity(2:end) - velocity(1:end-1))/dt;
disp = u(2:end-1);
external_force = (m*acceleration + c*v2 + k*disp)/(uforce*k);
load_err = max( external_force ) - min(external_force );
hi = load_err * time(2)^(-2); % ~1/20
