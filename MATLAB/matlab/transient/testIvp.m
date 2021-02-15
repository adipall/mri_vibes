%TESTIVP test solutionIvp
clear
[k,c,m] = getRandomModel();
uo = randn();
vo = randn();
n=100;
tf=1;
[u,velocity,acceleration]=solutionIvp(m,c,k,uo,vo,n,tf);
time = linspace(0,tf,n);
figure(1); plot(time,u);title('disp');
dt = time(2);
fdVelocity = (u(2:end) - u(1:end-1))/dt;
t2 = linspace(dt/2,tf-dt/2,n-1);
figure(2); plot(t2,fdVelocity,'k', time,velocity,'g');
title('velocity');
fdAcceleration = (velocity(2:end) - velocity(1:end-1))/dt;
figure(3); plot(t2,fdAcceleration,'k', time,acceleration,'g');
title('acceleration');
residual = m*acceleration+c*velocity+k*u;
assert( norm(residual,inf) < 1.e-12);
