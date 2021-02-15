%   function y = Ypm(KF,u,q,mv,nmax)
%   RADIAL MATHIEU FUNCTION OF THE SECOND KIND 
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients
%               -nmax= maximum order
%               -category:(even-even,even odd,odd even,odd odd)                                                         
%   OUTPUTS:    -y= vector of function values for all 'nmax' orders
%   'mv' is determined beforehand with function 'eig_Spm'
%   The Radial Mathieu Function is approximated by an expansion in
%   products of Bessel functions.
function y = Ypm(KF,u,q,mv,nmax)
v1=sqrt(q)*exp(-u);   
v2=sqrt(q)*exp(u);
y = zeros(nmax,1);
nCoeffs = size(mv,1);
vt = getVt(KF,nCoeffs);
if KF == 1  %even-even
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A0=Apm(1);
        yc = A0*besselj(0,v1)*bessely(0,v2);
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*besselj(jm,v1)*bessely(jm,v2);
        end
        r=fix(in/2);
        coef=((-1)^r)*sqrt(pi/2)/A0;
        y(k)=yc*coef;
    end
elseif KF == 2  %even-odd
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A1=Apm(1);
        yc = A1*(besselj(1,v1)*bessely(0,v2)+bessely(1,v2)*besselj(0,v1));
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*(besselj(jm,v1)*bessely(j,v2)+ ...
                      bessely(jm,v2)*besselj(j,v1));
        end 
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        y(k)=yc*coef;
    end
elseif KF == 3  %odd-even
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A2=Apm(1);
        yc = A2*(bessely(0,v2)*besselj(2,v1)-bessely(2,v2)*besselj(0,v1));    
        for j=2:nCoeffs         
            jm=fix(j-1);
            jp=fix(j+1);
            yc = yc + ((-1)^j)*Apm(j)*(bessely(jp,v2)*besselj(jm,v1)- ...
                      bessely(jm,v2)*besselj(jp,v1));
        end
        r=fix(in/2);
        coef=((-1)^r)*sqrt(pi/2)/A2;
        y(k)=yc*coef;
    end
elseif KF == 4  %odd-odd
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A1=Apm(1);
        yc = A1*(bessely(1,v2)*besselj(0,v1)-bessely(0,v2)*besselj(1,v1));
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*(bessely(j,v2)*besselj(jm,v1)- ...
                          bessely(jm,v2)*besselj(j,v1));
        end
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        y(k)=yc*coef;
    end
end
