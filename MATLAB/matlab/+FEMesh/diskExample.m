% cd code/Salinas/tools/matlab
n=8;
obj = FEMesh.disk();
[connectivity,coordinates] = obj.constructor(obj,n);
x = real(coordinates);
y = imag(coordinates);
coordinates = [x,y];
fig = 1;
obj.plotQuadMesh(connectivity,coordinates,fig);
boundary_nodes = find( sqrt( sum( coordinates.*coordinates,2)) > 1-1.e-8 );
assert( size(boundary_nodes,1) == 6*n );
elementType = 'Quad';
fexo = FEMesh.Exodus( elementType, connectivity, coordinates );
fexo.Title = 'Disk';
fexo.Filename = 'disk.g';
fexo.Time = 0;
fexo.Nodesets(1).ID = 1;
fexo.Nodesets(1).Nodes=boundary_nodes;
fexo.Nodesets(1).DistFactors=ones(6*n,1);%
fexo.Nodesets(1).Status=0;
fexo.Nodesets(1).Name='boundary';
filename = fexo.Filename;
cd loadExodusIntoMatlab
exo_put(filename,fexo);

