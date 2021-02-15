% 5. Find the worst constraints and corresponding node numbers
% in terms of rotational rigid body modes
p = 8;
casaConstraint,
tol = 7e4; % from visualizing r(:,4:6)
rotx = find(abs(r(:,4))>tol);
roty = find(abs(r(:,5))>tol); 
rotz = find(abs(r(:,6))>tol);
rowC= [rotx; roty; rotz];
% rowC=[10533,10535,17455,17471,17473,17475];
subC = constraint(rowC,:);
[row,col,val] = find(subC);
[~,perm] = sort(row);
col = col(perm);
n = size(col,1); % 20
subdomain = zeros(n,1);
for j=1:n,
    for i = 1:p,
        ordinal = find( globalRow{i} == col(j)-1 );
        if ~isempty(ordinal) 
            assert( subdomain(j) == 0 ); % dof on unique subdomains
            subdomain(j) = i;
        end
    end
end
nodes = zeros(n,1);
for j=1:n,
    i = subdomain(j);
    ordinal = find( globalRow{i} == col(j)-1 );
    nodes(j) = localMap{i}(ordinal,2);
end
%unique(nodes),
        %  s     s     m
nodes = [40656 40657 52257;...
         41330 41487 52417;...
         52665 52666 52450];