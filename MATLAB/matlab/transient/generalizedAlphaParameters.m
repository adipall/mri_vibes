%GENERALIZEDALPHAPARAMETERS function [alphaF,alphaM,beta,gamma]=
function [alphaF,alphaM,beta,gamma]=generalizedAlphaParameters(rho)
alphaF=rho/(rho+1);
alphaM=(2*rho-1)/(rho+1);
beta=0.25*(1-alphaM+alphaF)*(1-alphaM+alphaF);
gamma=0.5-alphaM+alphaF;
