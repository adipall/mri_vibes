% function [history,u,v,a,uh,vh,ah,time] = oneSdTest(id)
% function [history,u,v,a,uh,vh,ah,time] = oneSdTest(id)
[uo,vo,rho,mDamp,kDamp,~,~] = getSuiteParameters(id);
dtransient = sdIntegrator(mDamp,kDamp,rho);
levels = 21;
[~,deltaH]=getRhoH();
tf = (levels-1)*deltaH;
[K,M] = getKandM();%c=0 & uo=0 implies ao=0
history = dtransient.integrator(K,M,deltaH,levels,uo,vo);
time = linspace(0,tf,levels);
C = mDamp*M+kDamp*K;
[u,v,a]=solutionIvp(M,C,K,uo,vo,levels,tf);% exact
rampRate=0;
n=levels;
[uh,vh,ah]=generalizedAlpha(n,tf,M,C,K,rampRate,rho,uo,vo);% trap if rho =1
