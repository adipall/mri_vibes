% $ pwd 
% .../code/Salinas/tools/matlab
obj = FEMesh.icosahedron;
radius = sqrt( .5*(5+sqrt(5)));
coordinates = obj.coordinates/radius;
elementType = 'Tria3';
fexo = FEMesh.Exodus( elementType, obj.connectivity, coordinates );
fexo.Title = 'Icosahedron';
fexo.Filename = 'icosahedron.g';
fexo.Time = 0;
fexo.Sidesets(1).ID = 1;
numberSides = size(obj.connectivity,1);
inward = 2;
outward = 1;
% for side = 1:numberSides,
%    element = side;
%    ids = obj.connectivity(element,:);
%    centroid = sum( coordinates(ids,:))/3;
%    distrubtionFactor(side) = some_function( centroid );
% end
fexo.Sidesets(1).Elements=(1:numberSides);
fexo.Sidesets(1).Sides=ones(numberSides,1)*outward;
fexo.Sidesets(1).DistFactors=ones(numberSides,1);%
fexo.Sidesets(1).Status=0;
fexo.Sidesets(1).Name='exterior';
filename = fexo.Filename;
cd loadExodusIntoMatlab
%exo_put(filename,fexo);

