%   function joiningFactors = gpm(KF,parameter,mv,nmax)
%   mv= matrix of expansion coefficients from 'eig_Spm'
%   nmax= maximum order
%   category:(even even,even odd,odd even,odd odd)                                                                      
function joiningFactors = gpm(KF,parameter,mv,nmax)
joiningFactors= zeros(nmax,1);
point = pi/2;
nCoeffs = size(mv,1);
vt = getVt(KF,nCoeffs);
if KF == 1 || KF == 4
        Spm_pi2=Spm(1,point,mv,nmax);
else % KF 2,3
        dSpm_pi2=dSpm(2,point,mv,nmax);
end
if KF == 1 %even-even              
    for k=1:nmax
        A0=mv(1,k);     
        Spmk=Spm_pi2(k);                   
        r=fix(vt(k)/2);
        joiningFactors(k)=((-1)^r)*Spmk/(pi*A0); 
    end
elseif KF == 2  %even-odd         
    for k=1:nmax
        A1=mv(1,k);     
        dSpmk=dSpm_pi2(k); 
        r=fix((vt(k)-1)/2);
        joiningFactors(k)=-((-1)^r)*(1/(pi*sqrt(parameter)))*dSpmk/A1;
    end
elseif KF == 3  %odd-even             
    for k=1:nmax
        A2=mv(1,k);     
        dSpmk=dSpm_pi2(k);
        r=fix(vt(k)/2);
        joiningFactors(k)=((-1)^r)*(1/(parameter*pi))*dSpmk/A2;
    end
elseif KF == 4  %odd-odd
    for k=1:nmax
        A1=mv(1,k);     
        Spmk=Spm_pi2(k);
        r=fix((vt(k)-1)/2);
        joiningFactors(k)=((-1)^r)*(1/(sqrt(parameter)*pi))*Spmk/A1;   
    end
end
