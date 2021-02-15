string_parameters

mesh = gen_mesh(0, len, 101);
assert(verify_mesh(mesh));

mu_force = 0.75*len;
sigma_force = 0.05*len;
forcing = @(x)exp(-(x-mu_force).^2/(2*sigma_force^2)); % N / m
[K,M,F] = assemble(mesh, forcing);

force_time = @(t)(t==0);

np = size(mesh.p,1);

free = 2:(np-1);
n_free = np-2;

K_free = tension*K(free,free);
M_free = density*M(free,free);
F_free = F(free);

n_ew = 6;
[V,D] = eigs(K_free, M_free, n_ew, 'smallestreal');

%% transient
dt = .00001;
nt = 1000;

base_case.disp_vec = zeros(np,1);
base_case.force_time = force_time;
base_case.free = free;

direct = base_case;
direct.K = K_free;
direct.M = M_free;
direct.F = F_free;
direct.V = speye(n_free);
direct.label = 'Direct';

modal = base_case;
modal.K = D;
modal.M = speye(n_ew);
modal.F = F_free;
modal.V = V;
modal.label = 'Modal';

cases{1} = direct;
cases{2} = modal;

disp_vec = zeros(np,1);

plot_scale = [-1e-6 1e-6];
solve_second_order_system(mesh.p, dt, nt, plot_scale, cases, @(x,y)0, 0);
