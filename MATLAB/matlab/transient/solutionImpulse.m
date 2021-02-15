%SOLUTIONIMPULSE function u = solutionImpulse
% m u'' + c u' + k u =  k uforce,  t < T
%                    =  0          t > T
function u = solutionImpulse(m,c,k,uforce,T,num,tfinal)
if tfinal <= T,
   [u,~] = solutionConstantForce(m,c,k,uforce,num,tfinal);
   return;
end;
dt=tfinal/(num-1);  % Constant force phase
constantForceLevel=floor(T/dt)+1;%num = tfinal/dt + 1
constantForceTime=dt*(constantForceLevel-1);
[uc,~]=solutionConstantForce(m,c,k,uforce,constantForceLevel,constantForceTime);
[uT,vT]=solutionConstantForce(m,c,k,uforce,1,T); % cut off
d=.5*c/m;
omega=sqrt((k/m)-d*d);
bbar=(vT+uT*d)/(uT*omega);
t1=dt*constantForceLevel;% first discrete time after cut off
lag=t1-T;
t2=tfinal-T;% remaining time
n2=num-constantForceLevel;%!
time=linspace(lag,t2,n2);
sine=sin(omega*time);
decay=exp(-d*time);
cosine=cos(omega*time);
u2=uT*decay.*(cosine+sine*bbar);
u=[uc,u2];
