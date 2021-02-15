function mesh = gen_mesh(a, b, np)

mesh.p = linspace(a, b, np)';
mesh.e = [1 : (np-1); 2 : np]';
