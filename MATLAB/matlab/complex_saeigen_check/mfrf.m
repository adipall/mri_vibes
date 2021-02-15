% Modal FRF
function x=mfrf(freqs,forces)
% function x=mfrf(freqs,forces)
%
% compute modal FRF at frequencies in freqs, using matrices on disk
%  freqs is a list of frequencies (Hz) for the computation
%  forces has one column per frequency. Same number of rows as in matrices
%
% uses Kssr, Mssr, Cssrr, Cssri from Salinas

tmp=Mssr();
n=size(tmp,1);
m= tmp+tmp'-eye(n).*tmp;

k=Kssr();
k=k+k'-eye(size(k)).*k;

%cr = 0.000001*k;
%cr=zeros(size(k));
cr=Ccsrr();
%cr=cr+cr'-eye(size(cr)).*cr;

%ci=Ccsri();
%ci=ci+ci'-eye(size(ci)).*ci;

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
%       for FRF Ax-lam B x=F(t) 
B = [cr m;
          m Z];
 
A = [-k Z;
          Z m];
	  
 
% compute left, right eigenvectors from matlab
%[vr,lambda]=eig(full(A),full(B));
%[vl,lambda]=eig(full(A'),full(B'));

% Per Mike Ross's suggestion. vl and vr are left and right hand eigenvectors
[AA,BB,Q,Z,vr,vl]  = qz(full(A),full(B),'complex');

% check for diagonalization
ebtmp=vl';

%diag_matlab_A = vl'*A* vr;
%diag_matlab_B = vl'*B* vr;

%for i=1:n
%    ratio(i) = diag_matlab_A(i,i)/diag_matlab_B(i,i);
%end

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
% project force vector into state-space
    g = [forces(:,i);g_zero];
    
% do the superposition for this given frequency, w
% minus sign on alpha (know why) 
    for j=1:nmodes
      phi_dot_g = vl(:,j)'*g;
%	  phi_dot_g = dot(tmp(j,1:n),g(1:n));
%      if ( beta(j) ~= 0 && imag(alpha(j)/beta(j)) ~= 0 )
	    zj = phi_dot_g/(-alpha(j)+ii*w*beta(j));
	    w_out(:,i) = w_out(:,i) +  zj*vr(:,j);
%      end
    end
end    	
 
% extract solution "u" from "w"
x = w_out(1:n,:);

end
