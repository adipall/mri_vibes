string_parameters

mesh = gen_mesh(0, len, 101);
assert(verify_mesh(mesh));

forcing = @(x)(pi/len)^2*sin(pi*x/len); % N / m
[K,M,~] = assemble(mesh, forcing);

np = size(mesh.p,1);

free = 2:(np-1);

K = tension*K;
M = density*M;

%% eigenvalues

n_ew = 6;
[V,D] = eigs(K, M, n_ew, 'smallestreal');

figure(1), clf, hold on
for ie = 1:n_ew
    plot(mesh.p, V(:,ie))
    freq_i = sqrt(D(ie,ie))/(2*pi);
    fprintf('Mode %d frequency: %g (computed), %g (expected)\n', ...
        ie-1, freq_i, 220*(ie-1));    
end


