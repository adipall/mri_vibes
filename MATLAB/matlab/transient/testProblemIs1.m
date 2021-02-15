%test assumptions used in tests of software
% applied to model problem 1,  piston/piston1x1_ex2.inp
[K,M] = getKandM();
tol = 1.e-8;
assert( abs(K -1.5e+5) < tol );
assert( abs(M - .1) < tol );
[alpha,beta]=getProportionalDamping();
assert( abs(alpha - 446.635710098) < tol );
assert( abs(beta) < tol );
[rho,deltaH] = getRhoH();
assert( abs( rho - .9) < tol );
assert( abs( deltaH - 5e-4) < tol );