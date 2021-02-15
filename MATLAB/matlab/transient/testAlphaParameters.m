%TESTALPHAPARAMETERS test generalizedAlphaParameters
rho = rand();
assert( 0 <= rho );
assert( rho <= 1);
[alphaF,alphaM,beta,gamma]=generalizedAlphaParameters(rho);
e(1)= gamma - ( .5 - alphaM + alphaF );%order
assert( beta >= .25 + .5 * (alphaF-alphaM));
oneMmF = .5*(1-alphaM+alphaF);
e(2) = beta - oneMmF*oneMmF;
e(3) = 3*alphaF - alphaM - 1;
assert( norm(e,1) < 1.e-12 );
