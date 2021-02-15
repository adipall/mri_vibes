% function coords = sectorBoundary(rMin, rMax, thetaMin, thetaMax)
function coords = sectorBoundary(rMin, rMax, thetaMin, thetaMax)
azimuthal_number_nodes = 100;
dTheta = (thetaMax - thetaMin)/(azimuthal_number_nodes -1);
index = 1;
for axial = 1:azimuthal_number_nodes,
      theta = thetaMin + (axial-1)*dTheta;
      coords(index,1) = rMin * cos(theta);
      coords(index,2) = rMin * sin(theta);
      index = index + 1;
end
for axial = azimuthal_number_nodes:-1:1
      theta = thetaMin + (axial-1)*dTheta;
      coords(index,1) = rMax * cos(theta);
      coords(index,2) = rMax * sin(theta);
      index = index + 1;
end
coords(index,1) = rMin * cos(thetaMin);
coords(index,2) = rMin * sin(thetaMin);
