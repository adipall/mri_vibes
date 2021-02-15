% function y = dHpm1(category,u,q,mv,nmax)
% DERIVATIVE OF RADIAL MATHIEU FUNCTION OF THE THIRD KIND
% u    value of radial coordinate to compute function
% q    elliptical parameter (q > 0)
% mv   matrix of expansion coefficients from 'eig_Spm'
% nmax maximum order
% category (even even,even odd,odd even,odd odd)
function y = dHpm1(category,u,q,mv,nmax)
y1 = dJpm(category,u,q,mv,nmax);
y2 = dYpm(category,u,q,mv,nmax);
y = y1 + 1i*y2;
