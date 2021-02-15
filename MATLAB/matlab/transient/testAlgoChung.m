%TESTALGOCHUNG Chung&Hulbert's Generalized Alpha
clear
disp('contains idea that is probably wrong');
% Chung and Hulberts Generalized Alpha method
[k,c,m]=getRandomModel();
rampRate=randn();
rho= rand();% trapezoidal method?
tf=.01;
numRefine = 10;
err = zeros(numRefine,3);
h = zeros(numRefine,1); % h is the time step size
uo=0;
vo=0;
[alphaF,alphaM,~,~]=generalizedAlphaParameters(rho);
for i=1:numRefine,
    n=2^(i+1);
    [d,v,a]=generalizedAlpha(n,tf,m,c,k,rampRate,rho,uo,vo);
    time = linspace(0,tf,n);
    h(i) =  time(2);
    firstT = (1-alphaF)*h(i);  % time(1)=0, time(2)=h
    lastT = tf - alphaF*h(i);
    [u,velocity,~]=solutionRampedForce(m,c,k,rampRate,n,firstT,lastT);
    % todo: do not do this ...
    dd = (1-alphaF)*d(2:n) + alphaF * d(1:n-1);
    vv = (1-alphaF)*v(2:n) + alphaF * v(1:n-1);
    err(i,1)=norm(dd-u,inf);
    err(i,2)=norm(vv-velocity,inf);
    firstT = (1-alphaM)*h(i);
    lastT = tf - alphaM*h(i);
    [~,~,acceleration]=solutionRampedForce(m,c,k,rampRate,n,firstT,lastT);
    aa = (1-alphaM)*a(2:n) + alphaM * a(1:n-1);
    err(i,3)=norm(aa-acceleration,inf);
end
clf;
figure(1);
hold off;
scaleAccel = err(1,3)/h(1);
first = [h(1),h(numRefine)];

scaleSecond = sqrt( err(1,1)*err(1,2) )/(h(1)*h(1));
second = first.*first*scaleSecond;
range = first*scaleAccel;
loglog(h,err(:,1),'-.',h,err(:,2),'--',h,err(:,3),'-',first,range,'k-',first,second,'k-');
handle = legend('disp','velocity','accel','first order','second order','location','northwest');
handle.FontSize = 16;
% disp(diag(h.^(-1))*err);
