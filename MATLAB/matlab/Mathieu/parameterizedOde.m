function [t,y] = parameterizedOde(q,a,y0)
    tspan = [0 pi/2];
%    options = odeset('RelTol',1.e-11,'AbsTol',[1.e-10,1.e-10]);
    options = odeset('RelTol',1.e-4,'AbsTol',[1.e-4,1.e-4]);
    [t,y] = ode45(@nestedOde,tspan,y0,options);
    function dydt = nestedOde(t,y)
            dydt = [y(2); y(1)*( 2*q*cos(2*t)-a )];
    end
end
