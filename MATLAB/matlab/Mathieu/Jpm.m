%function values = Jpm(category,u,q,mv,nmax)
%RADIAL MATHIEU FUNCTION OF THE FIRST KIND   
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients from 'eig_Spm'
%               -nmax= maximum order
%               -category (even even, even odd, odd even,odd odd)                                                       
%    
%   The Radial Mathieu Function is approximated by an expansion of
%   products of Bessel functions. It is related to Spm(category,i*u,mv,nmax) by
%   Spm(category,i*u,mv,nmax)=sqrt(2*pi)*gpm(category,q,mv,nmax)*Jpm(category,u,q,mv,nmax)
function values = Jpm(category,u,q,mv,nmax)
values =[];
if category < 1 || category > 4
    return
end
if length(u) + length(q) > 2
    return
end
v1=sqrt(q)*exp(-u);   
v2=sqrt(q)*exp(u);
nCoeffs = size(mv,1);
vt = getVt(category,nCoeffs);
values = zeros(nmax,1);
if category == 1  %even-even
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A0=Apm(1);
        yc=A0*besselj(0,v1)*besselj(0,v2);
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*besselj(jm,v1)*besselj(jm,v2);
        end
        r=fix(in/2);
        coef=((-1)^r)*sqrt(pi/2)/A0;
        yc=yc*coef; 
        values(k)=yc;
    end
elseif category == 2 %even-odd
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A1=Apm(1);
        yc=A1*(besselj(1,v1)*besselj(0,v2) + besselj(1,v2)*besselj(0,v1)); 
        for j=2:nCoeffs
            jm=fix(j-1);
            yc = yc + ((-1)^jm)*Apm(j)*( besselj(jm,v1)*besselj(j,v2)+ ...
                      besselj(jm,v2)*besselj(j,v1) );
        end      
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        yc=yc*coef; 
        values(k)=yc;
    end
elseif category == 3 %odd-even
    for k=1:nmax  
        Apm=mv(:,k);  in=vt(k);
        A2=Apm(1);
        yc = A2*(besselj(0,v2)*besselj(2,v1)-besselj(2,v2)*besselj(0,v1));
        for j=2:nCoeffs
            jm=fix(j-1);
            jp=fix(j+1);
            yc = yc + ((-1)^j)*Apm(j)*(besselj(jp,v2)*besselj(jm,v1)- ...
                          besselj(jm,v2)*besselj(jp,v1));
        end
        r=fix(in/2);
        coef=((-1)^r)*sqrt(pi/2)/A2;
        yc=yc*coef;
        values(k)=yc;
    end
elseif category == 4 %odd-odd
    for k=1:nmax
        Apm=mv(:,k);
        in=vt(k);
        A1=Apm(1);
        yc = A1*(besselj(1,v2)*besselj(0,v1)-besselj(0,v2)*besselj(1,v1));
        for j=2:nCoeffs             
            jm=fix(j-1);    
            yc = yc + ((-1)^jm)*Apm(j)*(besselj(j,v2)*besselj(jm,v1)- ...
                          besselj(jm,v2)*besselj(j,v1));
        end
        r=fix((in-1)/2);
        coef=((-1)^r)*sqrt(pi/2)/A1;
        yc=yc*coef;
        values(k)=yc;
    end
end
