function x=cmfrf(freqs,forces)
% function x=cmfrf(freqs,forces)
% Complex Modal FRF
%
% compute modal FRF at frequencies in freqs, using matrices on disk
%  freqs is a list of frequencies (Hz) for the computation
%  forces has one column per frequency. Same number of rows as in matrices
%
% uses Kssr, Mssr

tmp=Mssr();
n=size(tmp,1);
m= tmp+tmp'-eye(n).*tmp;

k=Kssr();
k=k+k'-eye(size(k)).*k;

if exist('Cssrr.m','file')
   cr=Cssrr();
   % cr=cr+cr'-eye(size(cr)).*cr; % in acoustic we have C skew symetric
else
   cr=zeros(size(k));
   disp('Cssrr not found');
end

%if exist('Cssri.m','file')
%   ci=Cssri();
%   ci=ci+ci'-eye(size(ci)).*ci;
%else
%   ci=zeros(size(k));
%   disp('Cssri not found');
%end

if ( n ~= size(k,1) || n ~=  size(cr,1) ) 
  error('Matrix sizes don''t match.\n');
end

nf=max(size(freqs));

if size(freqs) ~= [n nf] 
  error('forces vector is wrong dimension.\n');
end

% form state-space equations

Z = zeros(length(k));
Z = sparse(Z);

% Form state-space
% NOTE: For eigenvalues, Ax = lam Bx
%       for FRF -Ax+i*lam B x=F(t) 

B = [cr m;
          m Z];
 
A = [-k Z;
          Z m];
	  
% compute left, right eigenvectors from matlab
% Per Mike Ross's suggestion. vl and vr are left and right hand eigenvectors

[AA,BB,Q,Z,vr,vl]  = qz(full(A),full(B),'complex');

% check for diagonalization
%diag_matlab_A = vl'*A* vr;
%diag_matlab_B = vl'*B* vr;
%ratio=diag(diag_matlab_A)./diag(diag_matlab_B);
%ratio=sort(ratio)/(2*pi);

% compute the alpha_i, beta_i
alpha=zeros(2*n,1);
beta=zeros(2*n,1);
for i=1:2*n
    alpha(i) = vl(:,i)' * A * vr(:,i);
    beta(i) = vl(:,i)' * B * vr(:,i);
end

% compute modal frf based on complex modal superposition

g_zero=zeros(n,1);
w_out=zeros(2*n,nf);
ii=sqrt(-1);
nmodes=n*2;

for i=1:nf
    f=freqs(i);
    w=2*pi*f;
    g = [forces(:,i);g_zero];
    for j=1:nmodes
      phi_dot_g = vl(:,j)'*g;
      zj = phi_dot_g/(-alpha(j)+ii*w*beta(j));
      w_out(:,i) = w_out(:,i) +  zj*vr(:,j);
    end
end    	
 
% extract solution "u" from "w"
x = w_out(1:n,:);

end
