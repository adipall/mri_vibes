function x=dfrf(freqs,forces)
% function x=dfrf(freqs,forces)
%
% compute direct FRF at frequencies in freqs, using matrices on disk
%  freqs is a list of frequencies (Hz) for the computation
%  forces has one column per frequency. Same number of rows as in matrices
%
% uses Kssr, Mssr, Cssrr, Cssri.

% tmp=Mssr();
tmp=mtilde();
n=size(tmp,1);
m= tmp+tmp'-eye(n).*tmp;

k=Kssr();
k=k+k'-eye(size(k)).*k;

% if exist('Cssrr.m','file')
if exist('ctilde.m','file')
   cr=ctilde();
   % cr=cr+cr'-eye(size(cr)).*cr; % in acoustic we have C skew symetric
else
   cr=zeros(size(k));
   % disp('Cssrr not found');
   disp('ctilde not found');
end

if exist('Cssri.m','file')
   ci=Cssri();
   ci=ci+ci'-eye(size(ci)).*ci;
else
   ci=zeros(size(k));
   disp('Cssri not found');
end

if ( n ~= size(k,1) || n ~=  size(cr,1) || n ~= size(ci,1)  ) 
  error('Matrix sizes don''t match.\n');
end

nf=max(size(freqs));

if size(freqs) ~= [n nf] then
  error('forces vector is wrong dimension.\n');
end

for i=1:nf
  f=freqs(i);
  w=2*pi*f;
  A=k - w^2*m + sqrt(-1)*w*(cr+sqrt(-1)*ci);
  x(:,i)=A\forces(:,i);
end


