% assumtions
% 1. Yes at each node the same constraint equation is used for all 3 dof.
% 2. Yes, the slave node is the first node in each constraint equation.
% simple_loft example only
function simpleloftFigure(x0,y0,z0,mpcOrder,mpcVector,newMpcOrder,newMpcVector)
coords = [x0,y0,z0];
master_nodes = mpcVector(2:mpcOrder(1),1);
plotQuad( master_nodes, coords,0); hold on;

numMpc = size(mpcOrder,1);
hi = [0;3*mpcOrder(1:3:numMpc-3)];
hi=cumsum(hi);
num_unique_mpc = numMpc/3;
slave_ids = hi + ones(num_unique_mpc,1);
slave_nodes = mpcVector(slave_ids,1);
flip = [1 2 4 3]; % simple_loft only, correction for orientation
plotQuad( slave_nodes(flip), coords,1);

projection = zeros(4,3);
id=0;
for i =1:num_unique_mpc,
    start = id+2;
    stop = id+mpcOrder(i);
    id = id + 3*mpcOrder(i);
    face = start:stop;
    sumCoords= sum( mpcVector(face,3) );
    affine_coordinates = mpcVector(face,3)'/sumCoords;
    master_nodes = mpcVector(face,1);
    projection(i,:) = affine_coordinates*coords(master_nodes,:);
end
plotQuad(flip,projection,2);
plot3( projection(:,1),projection(:,2),projection(:,3),'v')

p = zeros(4,3);
id = 0;
for i =1:4,
    n = newMpcOrder(i);
    start = id+2;
    stop = id+n;
    id=stop;
    face = start:stop;
    sumCoords= sum( newMpcVector(face,2) );
    scale = sign(sumCoords);
    coef = newMpcVector(face,2)'*scale;
    master_nodes = newMpcVector(face,1);
    p(i,:) = coef*coords(master_nodes,:);
end
disp(p),
plotQuad(flip,p,3);
plot3( p(:,1),p(:,2),p(:,3),'o')
hold off;
az=-176;el=-61;view(az,el);axis('off');