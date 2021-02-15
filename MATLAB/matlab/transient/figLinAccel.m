%LINACCELFIG accerations are lousy if rho < 1
clear
uo=0;vo=1;mDamp=0;kDamp=0; name = 'GeneralizedAlpha';rampRate=0;
rho = .5;
% mDamp=446.635710098;
[K,M] = getKandM();
dtransient = sdIntegrator(mDamp,kDamp,rho);
C=mDamp*M+kDamp*K; tf = 2*pi*sqrt(M/K);
numRefine = 10;
err = zeros(numRefine,6);
h = zeros(numRefine,1); % h is the time step size
[alphaF,alphaM,~,~]=generalizedAlphaParameters(rho);
for i=1:numRefine,
    n= 2^(i+1);
    [uh,vh,ah]=generalizedAlpha(n,tf,M,C,K,rampRate,rho,uo,vo);
    deltaH = tf/(n-1);
    history = dtransient.integrator(K,M,deltaH,n,uo,vo);
    h(i) =  deltaH;
    [u,v,a]=solutionIvp(M,C,K,uo,vo,n,tf);
    u=u';v=v';a=a';
    uh=uh';vh=vh';ah=ah';
    err(i,1)=norm(uh-u,inf);
    err(i,2)=norm(vh-v,inf);
    err(i,3)=norm(ah-a,inf);
    err(i,4)=norm(history(:,1)-u,inf);
    err(i,5)=norm(history(:,2)-v,inf);
    err(i,6)=norm(history(:,3)-a,inf);
end
loglog(h,err(:,1),'-.',h,err(:,2),'--',h,err(:,3),'-',h,err(:,4),'o',h,err(:,5),'+',h,err(:,6),'v');
handle = legend('disp','velocity','accel','SD disp','SD vel','SD accel','location','northwest');
title(name);
fontSize = 14;
handle.FontSize = fontSize;
xmin = h(numRefine);
xmax = h(1);
ymin = err(numRefine,1)*.99;
ymax = err(1,3);
axis([xmin,xmax,ymin,ymax]);
trix = [1.e-4,1.e-3,1.e-3,1.e-4];
triy = [1.e-7,1.e-7,1.e-5,1.e-7];
hold on;
loglog(trix,triy,'k');
one = text(3.1e-4, 1.6e-7,'1'); one.FontSize = fontSize;
two = text(8e-4, 1.e-6,'2'); two.FontSize = fontSize;
earx = [.3e-5,.3e-3,.3e-5,.3e-5];
eary = [4.,4.e2,4.e2,4.e0];
loglog(earx,eary,'k');
uno = text(.33e-5,.4e2,'1'); uno.FontSize = fontSize;
eine = text(.2e-4, 1.5e2,'1'); eine.FontSize = fontSize;
hold off;
% axis off;
% print -dpdf trape
