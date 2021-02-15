% function [V,E] = symmetricEig(J,KF)
function [V,E] = symmetricEig(J,KF)
if( KF == 1 ) % J(2,1) = 2 J(1,2)
    mean = sqrt(2)*J(1,2);
    T = J;
    T(2,1) = mean;
    T(1,2) = mean;
    [V,E] = eig(T); % T V = V Lambda
    V(1,:) = V(1,:)*sqrt(.5); % J D V = D V Lambda,
else
    [V,E] = eig(J);
end
