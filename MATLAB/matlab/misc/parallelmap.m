% function [localMap,globalRow]=parallelmap(name,numProc)
% example casaMap.m
% also see dofmap
function [localMap,globalRow]=parallelmap(name,numProc)
gids = loadGids(name,numProc);
fmap = loadMaps(numProc);
check = isVolumetic(fmap);
if sum(sum(check))>0,
    disp(check);
end
assert( sum(sum(check)) == 0);
localMap = cell(numProc,1);
for i = 1:numProc,
    localMap{i} =triplet(gids{i},fmap{i});
end
asetsize = zeros(numProc,1);
for i=1:numProc,
    asetsize(i) = size( localMap{i},1);
end
globalRow = cell(numProc,1);
allNodes = zeros(sum(asetsize),1);
first = 0;
for i=1:numProc,
    last = first + asetsize(i);
    allNodes(first+1:last) = localMap{i}(:,2);
    first = last;
end
nodeIds = unique(allNodes,'rows');

nnode = size(nodeIds,1);
activedof = zeros(nnode,3);
for i=1:numProc,
    m = asetsize(i);
    for j=1:m
        nodeNumber = find( localMap{i}(j,2) == nodeIds );
        k = localMap{i}(j,3);
        activedof(nodeNumber,k)=1;
    end
end
numdof = sum(sum(activedof));
nodeBegVec = zeros(nnode+1,1);
localDofVec=zeros(numdof,1);
nnz = 1;
for i=1:nnode,
   for j=1:3,
      if activedof(i,j)==1,
         localDofVec(nnz)=j;
         nnz = nnz+1;
      end
   end
   nodeBegVec(i+1) = nnz;
end
globalRow = cell(numProc,1);
for i=1:numProc,
    ordinal = 1;
    previous = -1;
    localdof = 0;
    while ordinal<=asetsize(i),
        nodeNumber = find( localMap{i}(ordinal,2) == nodeIds );
        if nodeNumber == previous
           localdof = localdof + 1;
        else
           nnz = nodeBegVec(nodeNumber);
           previous = nodeNumber;
           localdof = 0;
        end
        globalRow{i}(ordinal) = nnz+localdof;
        ordinal = ordinal + 1;
    end
end
