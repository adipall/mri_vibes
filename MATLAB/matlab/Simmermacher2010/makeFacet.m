rMIn = 1;  % Geometry:   a sector of an annulus
rMax = 2;  % inner and outer radii
thetaMin = 0;   % angles defining the sector
thetaMax = pi/4;

% nr is small and nTheta is big
nr = 4;
nTheta = 5;


[coords, connectivity ] = sectorMesh( rMIn, rMax, nr, thetaMin, thetaMax, nTheta );
nElement = (nr-1)*(nTheta-1);
figure(1);
color = 'r';
for element = 1:nElement,
   plotPoly( element, connectivity, coords, color);
end
%hold off;

% nr is big and nTheta is small
nr = 6;
nTheta = 2;



[coords, connectivity ] = sectorMesh( rMIn, rMax, nr, thetaMin, thetaMax, nTheta );
nElement = (nr-1)*(nTheta-1);
figure(1);
color = 'k';
for element = 1:nElement,
   plotPoly( element, connectivity, coords, color);
end
hold off;
axis('off');
