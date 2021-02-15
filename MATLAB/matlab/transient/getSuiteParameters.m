%function [uo,vo,rho,mDamp,kDamp,name,number] = getSuiteParameters(id)
function [uo,vo,rho,mDamp,kDamp,name,number] = getSuiteParameters(id)
% match tests to a list of integers
assert( size(id,1) == 3 && size(id,2) == 1);
% id(1) = ivp,  id(1)=0 homogeneous problem
% homogeneous ic is a trivial problem
assert( -1 < id(1) && id(1) < 3);
ivp = id(1);
% id(2) = numberical damping,  0: rho=1,  1: rho=.9
assert( -1 < id(2) && id(2) < 2);
setRho = id(2);
% id(3) = proportional damping,  0: alpha=0   1: alpha = 400
assert( -1 < id(3) && id(3) < 2);
setDamping = id(3);
switch ivp
    case 0
        uo = 0; vo = 0;
    case 1
        uo = 0; vo = 1;
    case 2
        uo = 1; vo = 0;
    otherwise
        error('ivp not in 0 1 2');
end
if setRho == 0,
    rho = 1;
    name = 'Trapezoidal';
else
    [rho,~] = getRhoH();
    name = 'Generalized\alpha';
end
if setDamping == 0,
    mDamp = 0; kDamp = 0;
else
    [mDamp,kDamp]=getProportionalDamping();
end
number = id(1) + 2*id(2) + 4*id(3);
