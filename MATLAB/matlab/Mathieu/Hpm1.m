% function values = Hpm1(category,u,q,mv,nmax)
% RADIAL MATHIEU FUNCTION OF THE THIRD KIND
% u    value of radial coordinate to compute function
% q    elliptical parameter (q > 0)
% mv   matrix of expansion coefficients from 'eig_Spm'
% nmax maximum order
% category (even-even,even-odd,odd-even,odd-odd)
function values = Hpm1(category,u,q,mv,nmax)
realPart = Jpm(category,u,q,mv,nmax);
imagPart = Ypm(category,u,q,mv,nmax);
values = realPart + 1i*imagPart;
