%   function y = dJpm(KF,u,q,mv,nmax)
%   DERIVATIVE OF RADIAL MATHIEU FUNCTION OF THE FIRST KIND
%   u  value of radial coordinate to compute function
%   q  elliptical parameter (q > 0)
%   mv matrix of expansion coefficients from 'eig_Spm'
%   nmax maximum order
%   category:(even-even, even odd, odd even, odd odd)
%   The Radial Mathieu function approximation is an expansion
%   of products of Bessel functions.
function y = dJpm(KF,u,q,mv,nmax)
v1=sqrt(q)*exp(-u);
v2=sqrt(q)*exp(u);
nCoeffs = size(mv,1);
vt = getVt(KF,nCoeffs);
y=zeros(nmax,1);
if KF == 1  % even even
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A0=Apm(1);
        yc = A0*(v1*besselj(1,v1)*besselj(0,v2)-v2*besselj(0,v1)*besselj(1,v2));
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*(v1*besselj(j,v1)*besselj(jm,v2)- ...
                v2*besselj(jm,v1)*besselj(j,v2));
        end
        r=fix(in/2);
        coef=((-1)^r)*sqrt(pi/2)/A0;
        y(k)=yc*coef;
    end
elseif KF == 2  %even odd
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A1=Apm(1);
        yc = A1*((v2-v1)*(besselj(0,v1)*besselj(0,v2)-besselj(1,v2)*besselj(1,v1))+...
            besselj(1,v1)*besselj(0,v2)-besselj(0,v1)*besselj(1,v2));
        for j=2:nCoeffs
            jm=fix(j-1);
            jm1=fix(2*j-1);
            yc = yc + ((-1)^jm)*Apm(j)*((v2-v1)*(besselj(jm,v1)*besselj(jm,v2)- ...
                besselj(j,v1)*besselj(j,v2))+jm1*(besselj(j,v1)*besselj(jm,v2)- ...
                besselj(jm,v1)*besselj(j,v2)));
        end
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        y(k)=yc*coef;
    end
elseif KF == 3  %odd-even
    for k=1:nmax
        Apm=mv(:,k);  in=vt(k);
        A2=Apm(1);
        yc = -A2*4*(besselj(0,v1)*besselj(0,v2)+ ...
            (cosh(2*u))*besselj(1,v1)*besselj(1,v2)-(1/v1)*besselj(1,v1)* ...
            besselj(0,v2)-(1/v2)*besselj(0,v1)*besselj(1,v2));
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^j)*Apm(j)*4*j*(besselj(jm,v1)*besselj(jm,v2)+ ...
                (cosh(2*u))*besselj(j,v1)*besselj(j,v2)-j*((1/v1)*besselj(j,v1)* ...
                besselj(jm,v2)+(1/v2)*besselj(jm,v1)*besselj(j,v2)));
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
        yc=A1*((v1+v2)*(besselj(0,v1)*besselj(0,v2)+besselj(1,v2)*besselj(1,v1))-...
            besselj(1,v1)*besselj(0,v2)-besselj(0,v1)*besselj(1,v2));
        for j=2:nCoeffs
            jm=fix(j-1);
            jm1=fix(2*j-1);
            yc = yc + ((-1)^jm)*Apm(j)*((v1+v2)*(besselj(jm,v1)*besselj(jm,v2)+ ...
                besselj(j,v2)*besselj(j,v1))-jm1*(besselj(j,v1)*besselj(jm,v2)+ ...
                besselj(jm,v1)*besselj(j,v2)));
        end
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        y(k)=yc*coef;
    end
end
