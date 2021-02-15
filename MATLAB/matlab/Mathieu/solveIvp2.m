parameter = 1;
eigenvalues = [-.4551386041074136, 4.371300982735086,  16.03383234035951];
n = 0;
if n==0,
    y0 = [1 0]; % even-even 
else % n=1 
    y0 = [0 1]; % even-odd? 
end
[t,y] = parameterizedOde(parameter,eigenvalues(1),y0);
plot(t,y(:,1),'k-',t,y(:,2),'b-.');
y(end,:),
