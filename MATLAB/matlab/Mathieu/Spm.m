%   function y = Spm(KF,v,mv,nmax)
%   ANGULAR MATHIEU FUNCTION 
%   INPUTS:     -v= value of angular coordinate in radians 
%               -mv= matrix of expansion coefficients from 'eig_Spm'
%               -nmax= maximum order 
%   category:(even even,even odd,odd even,odd odd)
function y = Spm(KF,v,mv,nmax)
nCoeffs = size(mv,1);
% vt = getVt(KF,nCoeffs);
y = zeros(nmax,1);
if KF == 1  %even-even
    for k=1:nmax
    Apm=mv(:,k); 
    yc = Apm(1);
    for j= 2:nCoeffs
        jp = fix(2*(j-1));
        yc= yc + Apm(j)*cos(jp*v);
    end
    y(k)=yc; 
    end
elseif KF == 2  %even-odd
    for k=1:nmax
        Apm=mv(:,k); 
        yc = 0; 
        for j = 1:nCoeffs
            yc = yc + Apm(j)*cos((2*j-1)*v);
        end
        y(k)=yc; 
    end
elseif KF == 3  %odd-even
    for k=1:nmax
        Apm=mv(:,k); 
        yc = 0;   
        for j = 1:nCoeffs
            yc = yc + Apm(j)*sin(2*j*v);
        end
        y(k)=yc;       
    end
elseif KF == 4  %odd-odd  
    for k=1:nmax
        Apm=mv(:,k);
        yc = 0;    
        for j = 1:nCoeffs             
            yc = yc + Apm(j)*sin((2*j-1)*v);
        end
        y(k)=yc; 
    end
end
