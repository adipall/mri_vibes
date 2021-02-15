% hypothesis: tripetTest passes, does not assume slave node is first
% in each mpc
clear
load('simple_loft.mat');
filename='simple_loft.txt';
[mpcVector,mpcOrder]=readMpc(filename);
threshold = .05;
numMpc = size(mpcOrder,1)/3;
residual = zeros(numMpc,1);
os=0;
m = nsteps;
coords = [x0,y0,z0] + [nvar01(:,m),nvar02(:,m),nvar03(:,m)];
% worst is 255
for k=1:numMpc,
    i= 3*(k-1)+1;
    n = mpcOrder(i);
    serialNode = zeros(n,1);
    for j = 1:n,
        serialNode(j) = find( node_num_map == mpcVector(os+j));
    end
    weight = mpcVector(os+1:os+n,3);
    R = [ones(n,1), coords(serialNode,:)];
    residual(k)= norm(R'*weight,1);
%     [U,S,V] = svd(R); % R' = V' S' U
%     s3(k) = S(3,3);
% In node face contact, R should have rank 3.
% But sometimes R has rank 4... when n > 5
%     s4(k) = S(4,4);
%     U = U(:,1:3);
%     w = weight - U*(U'*weight);
%     r2(k) = norm(R'*w,1);
    os = os + 3*n;
end
semilogy(residual,'o'),
title('|CR|');

normresidual=norm(residual,inf);
if  normresidual > threshold,
    message = ['residual norm =',num2str(normresidual)]
    disp(message);
end

% 1. visualize the mesh
% 2. How many nodes does each element
% have in common with a constraint?
% 3. What happens if nodes with no weight are filtered out?
% rank R          constraint
%   1              2 nonzeros  point
%    2            3 nonzeros line
%                 
