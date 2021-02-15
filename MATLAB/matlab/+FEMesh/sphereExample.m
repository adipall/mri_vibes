% example
% cd code/Salinas/tools/matlab
obj = FEMesh.cubedSphere();
%elementType = 'Tri3';
%fexo = FEMesh.Exodus( elementType, obj.connectivity, obj.coordinates);
state=1;
sphericalCoordinates = obj.getCoordinates(n);
[X,Y,Z] = obj.getCartesianCoordinates(sphericalCoordinates);
coords = [X,Y,Z];
conn  = obj.getConnectivity(n);
nf = size( conn, 1);
hist = zeros(nf,4); high = zeros(nf,1);  low = zeros(nf,1);
for f = 1:nf,
   ids = conn(f,1:4);
   vertices = coords( ids ,:);
   angles = 2;
   widths = obj.elementDiagnostic(vertices,angles);
   hist(f,:) = widths';
   high(f) = max(widths);
   low(f) = min(widths);
end
