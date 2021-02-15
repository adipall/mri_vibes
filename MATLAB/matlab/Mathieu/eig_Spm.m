% function [va,mv,vt]=eig_Spm(category,ellipticalParameter,numberCoefficients)
%  angular Mathieu functions
%  category (even-even,even-odd,odd-even,odd-odd)
%  va   vector of eigenvalues 'a', size [ncoeffs 1]
%  mv   square matrix of expansion coefficients
%  vt   indices size  [1 ncoeffs],
function [va,mv,vt]=eig_Spm(category,ellipticalParameter,numberCoefficients)
va = [];
mv = [];
vt = [];
if (category < 1) || (category > 4) || (ellipticalParameter <= 0)
    return
end
M = getMatrix(category,ellipticalParameter,numberCoefficients);
[u,d]   = symmetricEig(M,category);
% norms = modalResiduals(M,u,d);
% checkZero = max(norms),
[va,permutation]= sort(diag(d));% va is vector of eigenvalues 'a'
v= u(:,permutation);
vt = getVt(category,numberCoefficients);
% columnSumV=sum(v);
% mv = v*diag( ones(1,numberCoefficients)./columnSumV );
mv = v;
