%% function [xi,eta]=getElliptic(x,y,focalLength) Cartesian>Elliptic
function [xi,eta]=getElliptic(x,y,focalLength)
x = x/focalLength;
y = y/focalLength;
u=x.^2 + y.^2;     % = (cosh(2 xi) + cos(2 eta))/2
v=2*(x.^2 - y.^2); % = 1 + cosh(2 xi) cos(2 eta)
discriminant = u.*u - v + 1;
%if discriminant < 0,
   %discriminant = 0;
%end
cosh2xi =u+sqrt(discriminant);
twoXi = acosh( cosh2xi );
xi = .5*twoXi;
ch = cosh(xi);
sh = sinh(xi);
c = x./ch;
s = y./sh;
eta = atan2(s,c);
