string_parameters

mesh = gen_mesh(0, len, 101);
assert(verify_mesh(mesh));

forcing = @(x)0;
[K,M,~] = assemble(mesh, forcing);

K = tension*K;
M = density*M;

np = size(mesh.p,1);

free = 2:(np-1);
n_free = np-2;

K_free = K(free,free);
M_free = M(free,free);

%% prescribed acceleration (at ends)

disp_scale = 1;
disp_omega = 2*pi*330;
disp_time = @(t)disp_scale*sin(disp_omega*t);
disp_vec = zeros(np,1);
disp_vec(1) = 1;
disp_vec(end) = 1;

F_accel = disp_omega^2*M(free,:)*disp_vec - K(free,:)*disp_vec;

%% con mass

con_mass_factor = 1e4;
con_mass_size = con_mass_factor*len*density;

M_con = M;
M_con(1,1) = M_con(1,1) + con_mass_size;
M_con(end,end) = M_con(end,end) + con_mass_size;

n_ew = 8;
[Phi_con,D_con] = eigs(K, M_con, n_ew, 'smallestreal');

F_con = -disp_vec*disp_omega^2*con_mass_size;

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

con_mass.K = K;
con_mass.M = M_con;
con_mass.F = F_con;
con_mass.V = speye(np);
con_mass.force_time = disp_time;
con_mass.disp_vec = zeros(np,1);
con_mass.free = 1:np;
con_mass.v0_ind = [1 np];
con_mass.v0_val = disp_scale*disp_omega;
con_mass.label = 'Con Mass Direct';

con_mass_modal.K = D_con;
con_mass_modal.M = speye(n_ew);
con_mass_modal.F = F_con;
con_mass_modal.V = Phi_con;
con_mass_modal.force_time = disp_time;
con_mass_modal.disp_vec = zeros(np,1);
con_mass_modal.free = 1:np;
con_mass_modal.v0_ind = 1:n_ew;
con_mass_modal.v0_val = disp_scale*disp_omega*...
    Phi_con'*M_con*disp_vec;
con_mass_modal.label = 'Con Mass Modal';

cases{1} = direct;
cases{2} = con_mass;
cases{3} = con_mass_modal;

plot_scale = [-2 2];
solve_second_order_system(mesh.p, dt, nt, plot_scale, cases, @analytic_soln, 0);
