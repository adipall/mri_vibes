% function map=triplet(gid,fmap) [equation nodeIndex local_dof]
% See also dofmap
% example exTriplet.m
%   map(:,1) = 1:ndof always, map(:,[2,3]): origin
function map=triplet(gid,fmap)
dpn = getDofPerNode();
nnodes=size(gid,1);
ndofs=max(fmap);
assert( nnodes*dpn == size(fmap,1) );
map=zeros(ndofs,3);
dof=0;
for node=1:nnodes;
    base = (node-1)*dpn;
    nodeIndex = gid(node);
    for local_dof=1:dpn
        active = fmap(base+local_dof);
        if active>0
            dof = dof+1;
            map(dof,:)=[active nodeIndex local_dof];
        end
    end
end
