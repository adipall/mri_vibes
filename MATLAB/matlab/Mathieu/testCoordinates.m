%% test change of coordinates
clear
a = 2; 
b = 1;
focalLength = sqrt( a*a - b*b );
load major2minor1;  % coordinates p, and connectivity t
n = size(p,1);
[xi,eta] = getElliptic(p(:,1),p(:,2),focalLength);
x = zeros(n,2);
for i=1:n,
   [x(i,1),x(i,2)] = getCartesian(xi(i),eta(i),focalLength);
end
if norm(x-p,inf) < 1.e-12,
    disp('pass');
else
    disp('fail');
end
