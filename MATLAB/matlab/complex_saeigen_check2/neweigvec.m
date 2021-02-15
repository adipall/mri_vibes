function[uln,urn] = neweigvec(B,ul,ur,lam)
%% ------------------------------------------------------------------------
% function[uln,urn] = neweigvec(B,ul,ur,lam)
%
% This is a function to test for orthogonal eigenvectors, if there is a
% problem then we will orthogonalize.
%
% input: B (B matrix in linearization formulation)
%        ul (left eigenvector from dggev)
%        ur (right eigenvector from dggev)
%        lam (eigenvalues from dggev)
%
% output: uln (new left eigenvector)
%         urn (new right eigenvector)
%
% 12Nov08
%
% NOTE: not optimized for matlab. Instead, we plan to use this
% as a template for C++ code.
%% ------------------------------------------------------------------------

%% --- Parameters ---------------------------------------------------------
eps = 1.0e-6; % small number for comparisons
maxlam = max(abs(lam));
epslam = maxlam*eps;
%% ------------------------------------------------------------------------

%% --- algorithm ----------------------------------------------------------
Ble = length(B);

% Probably don't need to do this. We might be able to overwrite ul and ur,
% but for now.
uln = ul;
urn = ur;

for i = 1:Ble-1
    for j = (i+1):Ble
        if (abs(lam(i)-lam(j))<epslam)
            beta_11 = ul(:,i)'*B*ur(:,i);
            beta_12 = ul(:,i)'*B*ur(:,j);
            beta_21 = ul(:,j)'*B*ur(:,i);
            if(abs(beta_12)>abs(beta_11)*eps)
                urn(:,j) = ur(:,j)-beta_12/beta_11*ur(:,i);
                disp('Fix Right Vector');
                [i j]
            end
            %
            if(abs(beta_21)>abs(beta_11)*eps)
                uln(:,j) = ul(:,j) - beta_21/beta_11*ul(:,i); 
                disp('Fix Left Vector');
                [i j]
            end
        end
    end
end

% Do we know about the order from dggev. Can we assume that these come out
% in orthogonal pairs? If so, then we should not loop through everything.
