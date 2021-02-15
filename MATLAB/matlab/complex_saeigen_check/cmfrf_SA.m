function x=cmfrf_SA(freqs,forces,A,B,vl,vr)
% function x=cmfrf(freqs,forces)
% Complex Modal FRF
%
% compute modal FRF at frequencies in freqs, using matrices on disk
%  freqs is a list of frequencies (Hz) for the computation
%  forces has one column per frequency. Same number of rows as in matrices
%
% uses A, B, vl, vr

tmp=mtilde();
n=size(tmp,1);

nf=max(size(freqs));

if size(forces) ~= [n nf] 
  error('forces vector is wrong dimension.\n');
end

nA = size(vl,2);

% compute the alpha_i, beta_i
alpha=zeros(nA,1);
beta=zeros(nA,1);
for i=1:nA
    alpha(i) = vl(:,i)' * A * vr(:,i);
    beta(i) = vl(:,i)' * B * vr(:,i);
end

% compute modal frf based on complex modal superposition

g_zero=zeros(n,1);
w_out=zeros(2*n,nf);
ii=sqrt(-1);
nmodes=nA;

for i=1:nf
    f=freqs(i);
    w=2*pi*f;
    % g = [forces(:,i);g_zero];
    % You have to put the force vector in line with the  correct set of
    % equations.
    g = [g_zero;forces(:,i)];
    for j=1:nmodes
      phi_dot_g = vl(:,j)'*g;
      zj = phi_dot_g/(-alpha(j)+ii*w*beta(j));
      w_out(:,i) = w_out(:,i) +  zj*vr(:,j);
    end
end    	
 
% extract solution "u" from "w"
x = w_out(1:n,:);


% end
