function [count_fail,diag_fail] = diagonal_test(Gd,max_R);

%% ------------------------------------------------------------------
% (count_fail) = diagonal_test(Gd,Max_R)
%
% This is a script to test if a matrix is diagonal.
%
% input: Gd (Matrix to be tested)
%           Max_R (How close do you want the ratio of
%                       off diagonal/diagonal)
%
% output: count_fail (#of times it didn't pass)
%              diag_fail (Matrix with Values where there is off diagonal
%                               entries.
%
% 09Nov08, MRR
%% -------------------------------------------------------------------

%% --- Parameters ------------------------
eps = 1e-6; % small number
count_fail = 0.0; % Initialize value

%% --- Let er Rip
GdA = abs(Gd);
diag_fail=zeros(length(GdA),length(GdA));
count_fail = 0.0;
for k = 1:length(Gd)
    for l = 1:length(Gd)
        % Check if diagonal is near zero
        if (GdA(k,k) < eps)
            break
        end
        % Check for ratio
        ratio = GdA(k,l)/GdA(k,k);
        if (ratio > max_R && l ~= k)
            count_fail = count_fail + 1;
            diag_fail(k,l) = ratio;
        end
    end
end

%% --- May want to think about doing an eig on Gd and comparing diagonal
%% entries to the eigenvalues. 