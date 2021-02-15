%TESTSDINITIALSTEP regression test integrator level 1
accel = 0;
veloc = 1;
disp = 0;
[mDamp,kDamp]=getProportionalDamping();
[rho,deltaH] = getRhoH();
dtransient = sdIntegrator(mDamp,kDamp,rho);
[K,M] = getKandM();
rhs = dtransient.getRhs(K,M,deltaH,accel,veloc,disp);
rhsGold = 420.2331785504899;
assert( abs( rhs - rhsGold ) < 1e-10 );
dPlus=dtransient.linearSystem(K,M,deltaH,rhs);
dGold = 4.165918047660997e-4;
assert( abs(dPlus - dGold ) < 1.e-10);