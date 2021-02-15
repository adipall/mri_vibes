%% elliptic cylinder function
clear
      a = 2; b =1;  % geometry
load major2minor1;  % p, t
focalLength = sqrt( a*a - b*b );
[radii,angles] = getElliptic(p(:,1),p(:,2),focalLength);
s=struct('parameter',.1,'numberTerms',8,'category',1);
[axial,radial] = evaluateMathieu(s,angles,radii);
mode = 5;
height = axial(:,mode).*radial(:,mode);
trisurf(t,p(:,1),p(:,2),height);