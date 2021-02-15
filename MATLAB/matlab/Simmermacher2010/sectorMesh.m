% nR radial number of nodes 
% nTheta axial anumber of nodes 
% function [coords, connectivity ] = sectorMesh(rMin, rMax, nR, thetaMin, thetaMax, nTheta)
function [coords, connectivity ] = sectorMesh(rMin, rMax, nR, thetaMin, thetaMax, nTheta)
if nR <= 1,
   return
end
if nTheta <= 1,
   return
end
dr = (rMax - rMin)/(nR-1);
dTheta = (thetaMax - thetaMin)/(nTheta-1);
for radial = 1:nR
   r = rMin + (radial-1)*dr;
   for axial = 1:nTheta
      theta = thetaMin + (axial-1)*dTheta;
      index = axial + (radial-1) * nTheta;
      coords(index,1) = r * cos(theta);
      coords(index,2) = r * sin(theta);
   end
end
%  1     2     3     4    .... m = nTheta
% m+1   m+2   m+3    m+4       2m
% 2m+1  2m+2  2m+3             3m     
nQuad =  (nR-1) * (nTheta-1);
for radial = 1:nR-1,
   for axial = 1:nTheta-1
       quad = axial + (radial-1)*(nTheta-1);
       node = axial + (radial-1)*nTheta;
       connectivity(quad,:) = [ node,  node+nTheta, node+nTheta+1,node+1];
   end
end
