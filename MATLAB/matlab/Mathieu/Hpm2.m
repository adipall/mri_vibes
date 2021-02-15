%   function y = Hpm2(KF,u,q,mv,nmax)
%   RADIAL MATHIEU FUNCTION OF THE FOURTH KIND
%   INPUTS:     -u= value of radial coordinate to compute function
%               -q= elliptical parameter (q > 0)
%               -mv= matrix of expansion coefficients from 'eig_Spm'
%               -nmax= maximum order
%   category:(even even,even odd,odd even,odd odd)
function y = Hpm2(KF,u,q,mv,nmax)
y1 = Jpm(KF,u,q,mv,nmax);
y2 = Ypm(KF,u,q,mv,nmax);
y = y1 - 1i*y2;
