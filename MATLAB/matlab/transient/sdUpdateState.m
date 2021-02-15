%SDUPDATESTATE function [D_n,Vel_n,Accel_n]=
%  independent of model problem
function [D_n,Vel_n,Accel_n]=sdUpdateState(rho,deltaH,accel,veloc,disp,dPlus)
[~,~,newmarkBeta,newmarkGamma]=generalizedAlphaParameters(rho);
coeffOfAccelForVel=deltaH*((1-newmarkGamma)-newmarkGamma*(1-2*newmarkBeta)/(2*newmarkBeta));
coeffOfVelForVel=1-newmarkGamma/(newmarkBeta);
coeffOfDispForVel=newmarkGamma/(newmarkBeta*deltaH);
coefForVel=[coeffOfAccelForVel;coeffOfVelForVel;coeffOfDispForVel];
% mVelplus_temp=[mAccel_n,veloc,dPlus-disp]*coefForVel;

deltaSq = deltaH * deltaH;
coeffOfDispForAccel=1/(newmarkBeta*deltaSq);
coeffOfVelForAccel=-1/(newmarkBeta*deltaH);
coeffOfAccelForAccel=-(1-2*newmarkBeta)/(2*newmarkBeta);
coefForAccel=[coeffOfAccelForAccel;coeffOfVelForAccel;coeffOfDispForAccel];
Accel_n = [accel,veloc,dPlus-disp]*coefForAccel;
Velplus_temp=[accel,veloc,dPlus-disp]*coefForVel;
Vel_n = Velplus_temp;
D_n = dPlus;
