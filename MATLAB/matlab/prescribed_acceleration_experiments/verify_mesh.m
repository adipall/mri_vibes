function okay = verify_mesh(mesh)

okay = 1;

okay = okay && all(sort(mesh.p) == mesh.p);

ab = mesh.p(mesh.e);

okay = okay && all(ab(:,2) > ab(:,1));
