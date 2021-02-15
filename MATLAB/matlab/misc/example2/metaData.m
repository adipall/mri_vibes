% 1-based description of olio structure
olioBlock = cell(3 1);
olioBlock{1} = [1 2 3 5 6 8];
olioBlock{2} = [1 4 7 9];
olioBlock{3} = 3;
nodal_dof_per_block = [ 3 3 3 3 6 3 6 3 3];
numberSubdomains = 3;
nodesWithSixDof=[9 21 23 24];
o12 =[5 6 7 11   14   23 24 25 27 29];
o13 =      [11 12  15 16 17 18];
o23 =       11;
