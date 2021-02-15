%TESTSDSECONDSTEP regression test integrator level 2
accel = 0;
veloc = 1;
disp = 0;
[rho,deltaH] = getRhoH();
[mDamp,kDamp]=getProportionalDamping();
dtransient = sdIntegrator(mDamp,kDamp,rho);
[K,M] = getKandM();
rhs = dtransient.getRhs(K,M,deltaH,accel,veloc,disp);
rhsGold = 420.2331785504899;
assert( abs( rhs - rhsGold ) < .1 );
dPlus=dtransient.linearSystem(K,M,deltaH,rhs);
dGold = 4.165918047660997e-4;
[disp,veloc,accel] = dtransient.updateState(deltaH,accel,veloc,disp,dPlus);
rhs = dtransient.getRhs(K,M,deltaH,accel,veloc,disp);
rhsGold = 632.7400571931432;
assert( abs( rhs-rhsGold ) < 1.e-15 * rhsGold );