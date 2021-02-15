%EXAMPLERAMPED example exact solutions homogeneous
%  disp as well as ...
%  velocity, acceleration laid over their finite different approximations
clear
[m,c,k]=getRandomModel();
v = randn();
n = 100;
plus = n+1;
ti = 0;
tf = 1;
[u,~,~] = solutionRampedForce(m,c,k,v,plus,ti,tf);
time = linspace(0,tf,n);
figure(1); plot(time,u);title('disp');
dt = time(2);
velocity = (u(2:end) - u(1:end-1))/dt;
t2 = .5*(time(2:end) + time(1:end-1));
figure(2); plot(t2,velocity);title('velocity');
v2 = .5*(velocity(2:end) + velocity(1:end-1));
acceleration = (velocity(2:end) - velocity(1:end-1))/dt;
disp = u(2:end-1);
ramp = (m*acceleration + c*v2 + k*disp)/(v*k);
figure(3); plot(ramp); title('ramp');
