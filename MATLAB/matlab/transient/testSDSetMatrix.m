%TESTSDSETMATRIX regression test SD time integration
[alpha,beta]=getProportionalDamping();
[rho,deltaH] = getRhoH();
[K,M]=getKandM();
[alpha_f, alpha_m, newmarkBeta, newmarkGamma] = generalizedAlphaParameters(rho);
h2=deltaH*deltaH;
mScale=(1-alpha_m)/(newmarkBeta*h2)+alpha*(1-alpha_f)*newmarkGamma/(newmarkBeta*deltaH);
K_coeff = (1-alpha_f)+ beta*(1-alpha_f)*newmarkGamma/(newmarkBeta*deltaH);
C_coeff =(1-alpha_f)*newmarkGamma/(newmarkBeta*deltaH);
nscale = K_coeff/C_coeff;
nscaleGold = 2.506265664160402e-4;
assert( abs(nscale - nscaleGold) < 1.e-18 );
Cr = 0;  % This is weird.  I am not sure where Cr comes from.
conv = Cr + nscale*K;
conv=conv*C_coeff;
sigma=-1/mScale;
matrix = M - sigma * conv;
matrixGold = 0.10849084968820744;
assert( abs(matrix - matrixGold) < 1.e-15);
