function x=mfrf(freqs,forces)
% function x=mfrf(freqs,forces)
%
% Real Modal FRF
% compute modal FRF at frequencies in freqs, using matrices on disk
%  freqs is a list of frequencies (Hz) for the computation
%  forces has one column per frequency. Same number of rows as in matrices
%
% uses Kssr, Mssr
%
% Real solution with proportional damping

tmp=Mssr();
n=size(tmp,1);
m= tmp+tmp'-eye(n).*tmp;

k=Kssr();
k=k+k'-eye(size(k)).*k;

if ( n ~= size(k,1) )
  error('Matrix sizes don''t match.\n');
end

betak=1e-4;
nf=max(size(freqs));

if size(freqs) ~= [n nf] 
  error('forces vector is wrong dimension.\n');
end

[v,d]=eig(full(k),full(m));

% compute the alpha_i, beta_i
omegai=sqrt(diag(d));

% compute modal frf based on real modal superposition
g_zero=zeros(n,1);
w_out=zeros(n,nf);
ii=sqrt(-1);
nmodes=n;

for i=1:nf
    f=freqs(i);
    w=2*pi*f;
    g = [forces(:,i)];
    for j=1:nmodes
        phi_dot_g = v(:,j)'*g;
	zj = phi_dot_g/(omegai(j)*omegai(j) + ii*w*betak*omegai(j)*omegai(j) -w^2);
	w_out(:,i) = w_out(:,i) +  zj*v(:,j);
    end
end    	
 
x = w_out(1:n,:);

end
