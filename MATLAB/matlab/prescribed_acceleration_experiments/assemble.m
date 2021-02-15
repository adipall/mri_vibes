function [K,M,F] = assemble(mesh, forcing)

ne = size(mesh.e,1);
np = size(mesh.p,1);
n_entries = 4*ne;

K_entries = zeros(n_entries, 1);
M_entries = zeros(n_entries, 1);
rows = zeros(n_entries, 1);
cols = zeros(n_entries, 1);

F = zeros(np, 1);

for ie = 1:ne
    global_ind = 4*(ie-1) + (1:4);
    
    points = mesh.e(ie, :);
    loc_rows = points'*ones(1,2);
    loc_cols = ones(2,1)*points;

    a = mesh.p(points(1));
    b = mesh.p(points(2));
    
    locK = 1/(b-a)*[ 1 -1;
                    -1  1];
        
    locM = (b-a)/6*[2 1;
                    1 2];
    
    K_entries(global_ind) = locK(:);
    M_entries(global_ind) = locM(:);
    rows(global_ind) = loc_rows(:);
    cols(global_ind) = loc_cols(:);
    
    f1 = integral(@(x)(forcing(x).*(1-(x-a)./(b-a))), a, b);
    f2 = integral(@(x)(forcing(x).*((x-a)./(b-a))), a, b);
    F(points) = F(points) + [f1; f2];
end

K = sparse(rows, cols, K_entries, np, np);
M = sparse(rows, cols, M_entries, np, np);