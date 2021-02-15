%% ------------------------------------------------------------------------
% This script is to test the quadratic eigenvalue problem. 
%
% Here we will check that the diagonalization done in Salinas. In
% addition we will check if a reduced FRF will work.
%
% 24Oct08, MRR
%% ------------------------------------------------------------------------

clear all
close all

% This puts the files in this folder into your path
dr = pwd;
path(path,dr);

% You will want to change into the directory where you ran
% sa_eigen. In that folder there should be some .m files that are
% read into this code 
cd ../st_acoust_small.d

%% --- Parameters ---------------------------------------------------------
max_R = 1.0e-4;
nmodes_st = 32;
nmodes_ac = 32;
tnmodes = nmodes_st + nmodes_ac;
j = sqrt(-1);
eps = 10e-6;

%% ------------------------------------------------------------------------

%% --- Read in Values -----------------------------------------------------
alphai = Dggev_alphai;
alphar = Dggev_alphar;
beta = Dggev_beta;
A = Dggev_A';
B = Dggev_B';
% A = AMatrix';
% B = BMatrix';
vl = Dggev_vl';
vr = Dggev_vr';
% These are the onles after the orthogonalization in Salinas.
phi_l = phi_left';
phi_r = phi_right';


%% --- Create Eigenvectors and Eigenvalues. ------------------------------
for k = 1:length(alphar)
    if (beta(k) ~= 0)
        eigval(k,1) = (alphar(k) + j*alphai(k))/beta(k);
    else
        disp('what the beta?');
    end
end
eigdiag = sort(eigval); 

ul = dggev_evect(alphai,vl);
ur = dggev_evect(alphai,vr);

diagA = ul'*A*ur;
figure(1)
surf(abs(diagA));
diagB = ul'*B*ur;
l2 = diag(diagA)./diag(diagB);

[count_fail,diag_fail] = diagonal_test(diagA,max_R);

% Test if the eigenvectors of the rigid body modes are orthogonal
[eigdiag,ord] = sort(eigval); 
ul_vec=ul(:,ord);
ur_vec=ur(:,ord);
A_ord = A(ord,ord);
%diagA_v = ul_vec(:,13:length(A))'*A*ur_vec(:,13:length(A));     
figure(2)
spy(abs(diagA)>eps)
figure(3)
spy(abs(diagB)>eps)

% % Orthogonal Eigvenctors
[uln,urn] = neweigvec(B,ul,ur,eigval);

newA = uln'*A*urn;
newB = uln'*B*urn;
[count_fail_new,diag_fail_new] = diagonal_test(newA,max_R);

% Spy matrices for view of diagonalization
figure(4)
spy(abs(newA)>eps)

figure(5)
spy(abs(newB)>eps)

len = length(A);

% % ------ Test with new phi_left and phi_right from Salinas -----
zl = zeros(len,len);
zr=zeros(len,len);
for i = 1:len
    zl(:,i) = phi_l(:,2*i-1) + j.*phi_l(:,2*i);
    zr(:,i) = phi_r(:,2*i-1) + j.*phi_r(:,2*i);
end
        
salA = zl'*A*zr;
salB = zl'*B*zr;

[count_fail_sal,diag_fail_sal] = diagonal_test(salA,max_R);

% Spy matrices for view of diagonalization
figure(6)
spy(abs(salA)>eps)

figure(7) 
spy(abs(salB)>eps)

chkA = newA - salA;
chkB = newB - salB;

maxA = max(max(chkA));
if (maxA > eps)
    disp('ERROR  ERROR  ERROR')
    disp('Manoj: Honestly')
    %break
end


%% --- Play with FRF ----------------------------------------
% Reduction 
rn = 6:1:len; rn = rn'; % HERE YOU CAN PLAY WITH THE REDUCTION
newAr = uln(:,rn)'*A*urn(:,rn);
figure(10)
spy(abs(newAr)>eps)
newBr = uln(:,rn)'*B*urn(:,rn);
figure(11)
spy(abs(newBr)>eps)

% let us try the tilde FRF on the tilde system
freq = 10:1:100; freq = freq';
forces = zeros(tnmodes,length(freq));
% Apply a force at all nodes in all directions
forces(:,:) = 1;

x = dfrf_tilde(freq,forces);

figure(12)
semilogy(freq,abs(x(21,:)))

% let us try the modal FRF
xm=cmfrf_SA(freq,forces,A,B,uln(:,rn),urn(:,rn));

figure(13)
%semilogy(freq,abs(x(21,:)),freq,abs(xm(21,:)))
semilogy(freq,abs(x(21,:)))
hold on
semilogy(freq,abs(xm(21,:)),'--r');

% let us compare the ones in Salinas
xm_sal = cmfrf_SA(freq,forces,A,B,zl(:,rn),zr(:,rn));
% xm_sal = cmfrf_SA(freq,forces,A,B,phi_l(:,rn),phi_r(:,rn));

figure(14)
semilogy(freq,abs(xm_sal(21,:)),'--r');

figure(15)
semilogy(freq,abs(x(11,:)))
hold on
semilogy(freq,abs(xm_sal(11,:)),'--r');


%% Another Method for calculating the FRF
% 
% nf=max(size(freq));
% ii=sqrt(-1);
% g_zero=zeros(tnmodes,1);
% for i=1:nf
%     f=freq(i);
%     w=2*pi*f;
%     g = [g_zero;forces(:,i)];
%     phi_dot_g = uln(:,rn)'*g;
%     Hw = inv(-newAr+ii*w*newBr);
%     wo = Hw*phi_dot_g; % diagonalized space
%     xo = urn(:,rn)*wo; % Return to Reduced Space
%     %xm2(i) = xo(tnmodes+21);
%     xm2(i) = xo(21);
% end
% 
% figure(15)
% semilogy(freq,abs(x(21,:)),freq,abs(xm2))




%% --- Matlab Functions -----------------
% [V,D] = eig(A,B,'qz');
% [lamdiag,p] = sort(diag(D));
% Vec = V(:,p);
% [AA,BB,Q,Z,Vqz,Wqz] = qz(A,B);
% alphaqz = diag(AA);
% betaqz = diag(BB);
% lamQz = alphaqz./betaqz;
% 
% WqAVq = Wqz'*A*Vqz;
% figure(6)
% surf(abs(WqAVq))
% WqBVq = Wqz'*B*Vqz;
% 
% figure(7)
% spy(abs(WqAVq)>10e-6);
% [count_fail_WAV,diag_fail_WAV] = diagonal_test(WqAVq,max_R);

