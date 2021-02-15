function ev=dggev_ev(ai,ar,beta)
% function ev=dggev_ev(ai,ar,beta)
% convert the ai, ar, and beta from dggev to eigenvalues

n=max(size(ai));

for i=1:n
  ev(i)=(ar(i)+ai(i)*sqrt(-1))/beta(i);
end
