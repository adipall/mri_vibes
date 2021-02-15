%   functiony = dSpm(KF,v,mv,nmax)
%   DERIVATIVE OF ANGULAR MATHIEU FUNCTION
%   v= value of angular coordinate in radians
%   mv= matrix of expansion coefficients from 'eig_Spm'
%   nmax= maximum order
%   category (even even,even odd,odd even,odd odd)
function y = dSpm(KF,v,mv,nmax)
y = zeros(nmax,1);
nCoeffs = size(mv,1);
% vt = getVt(KF,nCoeffs);
if KF == 1  %even-even
    for k=1:nmax
        Apm=mv(:,k);
        yc = 0;
        for j = 2:nCoeffs
            jc = fix(2*(j-1));
            yc = yc - jc*Apm(j)*sin(jc*v);
        end
        y(k)=yc;
    end
elseif KF == 2  %even-odd
    for k=1:nmax
        Apm=mv(:,k);
        yc = 0;
        for j = 1:nCoeffs
            yc = yc - (2*j-1)*Apm(j)*sin((2*j-1)*v);
        end
        y(k)=yc;
    end
elseif KF == 3  %odd-even
    for k=1:nmax
        Apm=mv(:,k);
        yc = 0;
        for j = 1:nCoeffs
            yc = yc + 2*j*Apm(j)*cos(2*j*v);
        end
        y(k)=yc;
    end
elseif KF == 4  %odd-odd
    for k=1:nmax
        Apm=mv(:,k);
        yc = 0;
        for j = 1:nCoeffs
            yc = yc + (2*j-1)*Apm(j)*cos((2*j-1)*v);
        end
        y(k)=yc;
    end
end
