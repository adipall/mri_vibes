rMin = 1.5;  % Geometry:   a sector of an annulus
rMax = 2;  % inner and outer radii
thetaMin = 0;   % angles defining the sector
thetaMax = pi/4;
whatITSsees = sectorBoundary( rMin, rMax, thetaMin, thetaMax);

theta = .5*(thetaMin + thetaMax);
radius = .5*(rMin + rMax );
center = [cos(theta), sin(theta) ]* radius;

figure(1);
color = 'r';
plot( whatITSsees(:,1), whatITSsees(:,2), 'b');
hold on;
plot( center(1),center(2),'bs');
whatExodusSees = quad( rMin, rMax, thetaMin, thetaMax);
plot( whatExodusSees(:,1), whatExodusSees(:,2), 'r');
axis([.5, 2, 0, 1.5]);
axis('square');

centroid = sum( whatExodusSees(1:4,:))/4;
plot( centroid(1), centroid(2), 'ro');

hold off;
axis('off');
