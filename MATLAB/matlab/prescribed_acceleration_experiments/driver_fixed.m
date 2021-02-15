string_parameters

mesh = gen_mesh(0, len, 101);
assert(verify_mesh(mesh));

forcing = @(x)(pi/len)^2*sin(pi*x/len); % N / m
[K,M,F] = assemble(mesh, forcing);

np = size(mesh.p,1);

free = 2:(np-1);

K_free = tension*K(free,free);
M_free = density*M(free,free);
F_free = F(free);

%% static load
u_free = K_free\F_free;
u = zeros(np,1);
u(free) = u_free;
u_exact = sin(pi*mesh.p/len)/tension;

figure(1), clf
plot(mesh.p, u)
fprintf('Static load relative error: %g\n', norm(u-u_exact)/norm(u_exact))

%% eigenvalues

n_ew = 6;
[V,D] = eigs(K_free, M_free, n_ew, 'smallestreal');

figure(2), clf, hold on
for ie = 1:n_ew
    mode_i = zeros(np,1);
    mode_i(free) = V(:,ie);
    plot(mesh.p, mode_i)
    freq_i = sqrt(D(ie,ie))/(2*pi);
    fprintf('Mode %d frequency: %g (computed), %g (expected)\n', ...
        ie, freq_i, 220*ie);    
end


