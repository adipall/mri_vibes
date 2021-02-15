% example Generalized Alpha method
clear
[k,c,m]=getRandomModel();
rampRate=randn();
rho=rand();
tf=1;
uo=0;
vo=0;
refinements = 10;
err=zeros(refinements,3);
for i=1:refinements,
    n=2^i;
    [d,v,a]=generalizedAlpha(n,tf,m,c,k,rampRate,rho,uo,vo);
    time = linspace(0,tf,n);
    h =  time(2);
    resid = m*a+c*v+k*d-time*rampRate*k;
    quality = norm(resid,1);
    ti =0;
    plus = n+1;
    [u,velocity,acceleration] = solutionRampedForce(m,c,k,rampRate,plus,ti,tf);%Ag
    err(i,1)=norm(d-u,inf);
    err(i,2)=norm(v-velocity,inf);
    err(i,3)=norm(a-acceleration,inf);
end
