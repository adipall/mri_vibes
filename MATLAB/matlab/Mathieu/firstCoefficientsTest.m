% test Fourier coefficients of angular Mathieu functions
clear
n = 4;
nmax = 7;
[coefficients,q] = first4Coefficients();
computed = zeros(16,1);
last = 0;
for category = 1:4, 
    [~,mv,~]=eig_Spm(category,q,nmax);
    first = last + 1;
    last = last + 4;
    computed(first:last) =  mv(1:4,1);
end
if norm(computed - coefficients,inf) > 1.e-8,
    disp('Fail');
else
    disp('Pass');
end