% test Fourier coefficients of angular Mathieu functions
clear
category = 1;            % Gutierrez Vega Example 2
q=.1;                    % elliptical parameter
n = 6;                   % 0:2:10,  1-based 
coefficients = coefficientCe10();
numberCoefficients = size(coefficients,1); % 14
step = 0;
err = ones(1,9);
for nmax = n:numberCoefficients,
    [va,mv,vt]=eig_Spm(category,q,nmax);
    f = mv(:,n);
    if mv(n,n) < 0,
        f = -f;
    end
    f = f/( sqrt(2) * norm(f) );
    step = step + 1;
    err(step) = norm( f - coefficients(1:nmax) );
end
bound = [ 1.e-5,1e-8,1.e-10, ones(1,6)*(1.e-13)];
if( max(err/bound) > 1.)
     disp('Fail');
else
     disp('Pass');
end
