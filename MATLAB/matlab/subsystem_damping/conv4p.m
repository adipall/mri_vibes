function [fs,beta,kt]=conv4p(x,R,S,phimax)
%function [fs,beta,kt]=conv4p(x,R,S,phimax)
% converts from the old 4 parameter model to a new.
% based on Journal of Applied Mechanics, Sept 2005, vol 72. p. 755

fs=R*phimax^(x+2)/(x+2) + S*phimax;

beta=S/(R*phimax^(x+1)/(x+1));

c=beta+(x+1)/(x+2);

kt=fs*(1+beta)/(phimax*c);

fprintf('fs=%g, beta=%g, kt=%g\n',fs,beta,kt);
