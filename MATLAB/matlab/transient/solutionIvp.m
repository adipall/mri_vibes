%SOLUTIONIVP function [u,v,a] =
%ma+cv+k=0, u(0)=d_o, v(0)=v_o
function [u,v,a]=solutionIvp(m,c,k,uo,vo,n,tf)
d=.5*c/m;
omega=sqrt((k/m)-d*d);
time = linspace(0,tf,n);
decay  = exp( -d*time );
cosine = cos( omega * time );
sine   = sin( omega * time );
sigmau =(vo+d*uo)/omega;
sigmav =-(k*uo/m+d*vo)/omega;
ao=omega*sigmav-d*vo;
sigmaa=-(k*vo/m+d*ao)/omega;
u = decay.*(cosine*uo+sine*sigmau);
v = decay.*(cosine*vo+sine*sigmav);
a = decay.*(cosine*ao+sine*sigmaa);
