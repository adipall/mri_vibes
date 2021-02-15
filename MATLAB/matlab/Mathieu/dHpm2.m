%   function y = dHpm2(KF,u,q,mv,nmax)
%   DERIVATIVE OF RADIAL MATHIEU FUNCTION OF THE FOURTH KIND    
%   INPUTS:     -u= value of radial coordinate to compute function 
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients from 'eig_Spm'
%               -nmax= maximum order 
%               -KF= function code:  KF=1 even-even
%                                    KF=2 even-odd
%                                    KF=3 odd-even
%                                    KF=4 odd-odd      
%   OUTPUTS:    -y= vector of derivative values for all 'nmax' orders 
function y = dHpm2(KF,u,q,mv,nmax)
y1 = dJpm(KF,u,q,mv,nmax);
y2 = dYpm(KF,u,q,mv,nmax);
y = y1 - 1i*y2;

