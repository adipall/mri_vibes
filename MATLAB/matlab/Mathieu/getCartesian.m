% function [x,y] = getCartesian(xi,eta,focalLength) Elliptic -> Cartesian
function [x,y] = getCartesian(xi,eta,focalLength)
x = [];
y = [];
if eta < -pi || eta > pi,
  return;
end
x = focalLength*cosh(xi)*cos(eta);
y = focalLength*sinh(xi)*sin(eta);
