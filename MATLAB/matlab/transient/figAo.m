%AOFIG initial acceleration
%This is the figure about the initial acceleration
ivp = 2;
setRho = 0;
damping = 0;
id = [ivp,setRho,damping]';
[history,u,v,a,uh,vh,ah,time] = oneSdTest(id);
[~,~,~,~,~,name,number] = getSuiteParameters(id);
plot(time,history(:,3),time,a,time,ah);
xlabel('time'); ylabel('accel');
handle=legend('SD','Exact',name,'location','north');
handle.FontSize = 16;
