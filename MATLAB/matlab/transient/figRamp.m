%RAMPFIG time dependent forces are lousy
clear
uo=0;vo=0;rampRate=1;kDamp=0;mDamp=0;% mDamp=446.635710098;
rho = 1; % .5;
name = 'Trapezoidal'; % 'GeneralizedAlpha';
[K,M] = getKandM();
dtransient = sdIntegrator(mDamp,kDamp,rho);
dtransient.rampRate = rampRate;
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
    ti=0;levels=n+1;
    [u,v,a]=solutionRampedForce(M,C,K,rampRate,levels,ti,tf);
    u=u';v=v';a=a';
    uh=uh';vh=vh';ah=ah';
    err(i,1)=norm(uh-u,inf);
    err(i,2)=norm(vh-v,inf);
    err(i,3)=norm(ah-a,inf);
    err(i,4)=norm(history(:,1)-u,inf);
    err(i,5)=norm(history(:,2)-v,inf);
    err(i,6)=norm(history(:,3)-a,inf);
end
figure(1);
loglog(h,err(:,1),'-.',h,err(:,2),'--',h,err(:,3),'-',h,err(:,4),'+-.',h,err(:,5),'+--',h,err(:,6),'-+');
handle = legend('disp','velocity','accel','SD disp','SD vel','SD accel','location','westoutside');
title(name);
fontSize = 14;
handle.FontSize = fontSize;
xmin = h(numRefine);
xmax = h(1);
ymin = err(numRefine,1)*.99;
ymax = err(1,6);
axis([xmin,xmax,ymin,ymax]);
trix = [1.e-4,1.e-3,1.e-3,1.e-4];
triy = [1.e-7,1.e-7,1.e-5,1.e-7];
hold on;
loglog(trix,triy,'k');
one = text(3.1e-4, 1.6e-7,'1'); one.FontSize = fontSize;
two = text(8e-4, 1.e-6,'2'); two.FontSize = fontSize;
hold off;
% axis off;
% print -dpdf ramp
