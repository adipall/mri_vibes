% test that mpcs come in triplets,
% 3 dof per node,
% using the same mpc for each dof
clear
filename='simple_loft.txt';
[mpcVector,mpcOrder]=readMpc(filename);
% Do the mpcs come in 3's?
threshold = 1.e-12;
numMpc = size(mpcOrder,1);
xorder = mpcOrder(1:3:end);
yorder = mpcOrder(2:3:end);
assert( norm( yorder - xorder,1)== 0 ),
clear yorder
zorder = mpcOrder(3:3:end);
assert( norm( zorder - xorder,1)== 0 ),
clear xorder zorder
os=0;
for i=1:3:numMpc,
    n = mpcOrder(i);
    node = zeros(n,3);
    weight = zeros(n,3);
    for j=1:3,
        node(:,j)   = mpcVector(os+1+(j-1)*n:os+j*n,1);
        weight(:,j) = mpcVector(os+1+(j-1)*n:os+j*n,3);
        if j>1,
            assert( norm( node(:,j) - node(:,1),1)==0);
            assert( norm( weight(:,j) - weight(:,1),1)< threshold);
        end
    end
    os = os + 3*n;
    clear node weight
end
subplot(2,1,1);plot(mpcOrder);
subplot(2,1,2);plot(mpcVector(:,3),'.')
disp('Yes at each node the same constraint equation is used for all 3 dof');
conclusion=testMasterSlaveOrder(mpcOrder,mpcVector);
disp(conclusion);

