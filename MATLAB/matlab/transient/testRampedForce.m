%TESTRAMPEDFORCE test solutionRampedForce
clear
[k,c,m] = getRandomModel();
rampRate= randn();
n=100;
tf=1;
ti=0;
gaa = n + 1;
[u,velocity,acceleration] = solutionRampedForce(m,c,k,rampRate,gaa,ti,tf);
time = linspace(0,tf,n);
figure(1);
subplot(3,1,1);plot(time,u); ylabel('disp');
dt = time(2);
fdVelocity = (u(2:end) - u(1:end-1))/dt;
t2 = linspace(dt/2,tf-dt/2,n-1);
subplot(3,1,2); plot(t2,fdVelocity,'k', time,velocity,'g');
ylabel('veloc');
fdAcceleration = (velocity(2:end) - velocity(1:end-1))/dt;
subplot(3,1,3); plot(t2,fdAcceleration,'k', time,acceleration,'g');
ylabel('accel');
residual = m*acceleration+c*velocity+k*u-time*rampRate*k;
assert( norm(residual,inf) < 1.e-12);
