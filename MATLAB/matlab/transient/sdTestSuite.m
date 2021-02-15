%SDTESTSUITE test initial value problems, f=0
for ivp = 1:2,
for setRho = 0:1,
for damping = 0:1,
    id = [ivp,setRho,damping]';
    [history,u,v,a,uh,vh,ah,time] = oneSdTest(id);
    [~,~,~,~,~,name,number] = getSuiteParameters(id);
    figure(number);
    subplot(5,1,1);plot(time,history(:,1),time,u,time,uh);
    ylabel('disp');
    subplot(5,1,2);plot(time,history(:,2),time,v,time,vh);
    ylabel('veloc');
    subplot(5,1,3);plot(time,history(:,3),time,a,time,ah);
    ylabel('accel');
    u =u';  v=v'; uh = uh'; vh = vh';
    figure(number);
    subplot(5,1,4);plot(time,history(:,1)-u,time,uh-u);
    ylabel('disp');
    subplot(5,1,5);plot(time,history(:,2)-v,time,vh-v);
    ylabel('veloc');
    handle=legend('SD','Exact',name,'Disp Error','Vel Error','location','westoutside');
    handle.FontSize = 16;
end
end
end
