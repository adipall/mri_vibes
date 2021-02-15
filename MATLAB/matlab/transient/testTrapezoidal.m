%TESTTRAPEZOIDAL test generalizedA(rho=1)
clear
uo=0;vo=1;rho=1;mDamp=0;kDamp=0; name = 'Trapezoidal';rampRate=0;
[K,M] = getKandM();
C=mDamp*M+kDamp*K; tf = 2*pi*sqrt(M/K);
numRefine = 10;
err = zeros(numRefine,3);
h = zeros(numRefine,1); % h is the time step size
[alphaF,alphaM,~,~]=generalizedAlphaParameters(rho);
for i=1:numRefine,
    n= 2^(i+1);
    [uh,vh,ah]=generalizedAlpha(n,tf,M,C,K,rampRate,rho,uo,vo);
    time = linspace(0,tf,n);
    h(i) =  time(2);
    [u,v,a]=solutionIvp(M,C,K,uo,vo,n,tf);
    err(i,1)=norm(uh-u,inf);
    err(i,2)=norm(vh-v,inf);
    err(i,3)=norm(ah-a,inf);
end
loglog(h,err(:,1),'-.',h,err(:,2),'--',h,err(:,3),'-');
handle = legend('disp','velocity','accel','location','northwest');
title(name);
handle.FontSize = 16;
