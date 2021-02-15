string_parameters

mesh = gen_mesh(0, len, 101);
assert(verify_mesh(mesh));

forcing = @(x)0;
[K,M,~] = assemble(mesh, forcing);

np = size(mesh.p,1);

free = 2:(np-1);
n_free = np-2;

K_free = tension*K(free,free);
M_free = density*M(free,free);

n_ew = 6;
[Phi,D] = eigs(K_free, M_free, n_ew, 'smallestreal');

%% prescribed acceleration (at ends)

disp_omega = 2*pi*330;
disp_time = @(t)sin(disp_omega*t);
disp_vec = zeros(np,1);
disp_vec(1) = 1;
disp_vec(end) = 1;

F_accel = disp_omega^2*density*M(free,:)*disp_vec - ...
        tension*K(free,:)*disp_vec;

%% small eigenvalue problem for augmented matrix

Psi = ones(n_free, 1);

P = [Phi, Psi];
M_hat = P'*M_free*P;
K_hat = P'*K_free*P;

% symmetrize so that eig will return mass-normalized V
M_hat = 0.5*(M_hat + M_hat');
K_hat = 0.5*(K_hat + K_hat');

[V,Theta] = eig(K_hat, M_hat);

n_augment = size(P,2);

%% transient
dt = .00001;
nt = 1000;

base_case.disp_vec = disp_vec;
base_case.force_time = disp_time;
base_case.free = free;

direct = base_case;
direct.K = K_free;
direct.M = M_free;
direct.F = F_accel;
direct.V = speye(n_free);
direct.label = 'Direct';

modal = base_case;
modal.K = D;
modal.M = speye(n_ew);
modal.F = F_accel;
modal.V = Phi;
modal.label = 'Modal';

augment = base_case;
augment.K = Theta;
augment.M = speye(n_augment);
augment.F = F_accel;
augment.V = P*V;
augment.label = 'Augmented';

cases{1} = direct;
cases{2} = modal;
cases{3} = augment;

plot_scale = [-2 2];
solve_second_order_system(mesh.p, dt, nt, plot_scale, cases, @analytic_soln, 0);
