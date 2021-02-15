%SOLUTIONRAMPEDFORCE function [u,v,a]=
% mu''+cu'+ku=vkt
% The parameter v specifies the force, & has units of velocity.
function [u,velocity,acceleration]= solutionRampedForce(m,c,k,v,n,ti,tf)
dd=.5*c/m;
omega=sqrt((k/m)-dd*dd);
beta=c/k;
assert(beta*sqrt(k/m)<2);
gamma = v*beta;
sigma = (gamma*dd-v)/omega;
time = linspace(ti,tf,n-1);
decay  = exp(-dd*time);
cosine = cos(omega*time);
sine   = sin(omega*time);
eine = ones(1,n-1);
u = decay.*(cosine*gamma+sine*sigma)+(time-eine*beta)*v;
s = (beta*k-m*dd)/(m*omega);
velocity = v*(eine-decay.*(cosine+sine*s));
cc = dd - omega*s;
ss = dd*s + omega;
acceleration= v*(decay.*(cosine*cc+sine*ss));

