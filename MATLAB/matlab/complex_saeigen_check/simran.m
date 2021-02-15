%% ------------------------------------------------------------------------
% This script is to test the quadratic eigenvalue problem. 
%
% 24Oct08, MRR
%% ------------------------------------------------------------------------

close all

%% --- Parameters ---------------------------------------------------------
max_R = 1.0e-4;
%% ------------------------------------------------------------------------

%% --- Read in Values -----------------------------------------------------
alphai = Dggev_alphai;
alphar = Dggev_alphar;
beta = Dggev_beta;
A = Dggev_A';
B = Dggev_B';
vl = Dggev_vl';
vr = Dggev_vr';


%% --- Create Eigenvectors and Eigenvalues.
for k = 1:length(alphar)
    if (beta(k) ~= 0)
        eigval(k,1) = (alphar(k) + i*alphai(k))/beta(k);
    else
        disp('what the beta?');
    end
end
eigdiag = sort(eigval); 

ul = dggev_evect(alphai,vl);
ur = dggev_evect(alphai,vr);

[uln,urn] = neweigvec(A,ul,ur,eigval);

newA = uln'*A*urn;
newB = uln'*B*urn;
