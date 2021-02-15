clear
numberSubdomains = 4;
subdAd = getSubdAd();
adjacentLocalId = getAdjacentLocalId();
[matchListSize,matchList,owner]= getAdjacency(numberSubdomains);
[forceVector, nodalSum, displacementVector] = getFeVectors(numberSubdomains);
for subdomain = 1:numberSubdomains,
    myNodalSum = forceVector{subdomain};
    numberNeighbors = size(matchListSize{subdomain},1) - 1;
    for i = 1:numberNeighbors,
        neighbor = subdAd{subdomain}(i);
        start = matchListSize{subdomain}(i) + 1;
        stop = matchListSize{subdomain}(i+1);
        numberShared = stop - start + 1;
        eine = int32( ones(numberShared,1) );
        sharedIds = matchList{subdomain}(start:stop) + eine;
        j = adjacentLocalId{subdomain}(i);
        first = matchListSize{neighbor}(j) + 1;
        last  = matchListSize{neighbor}(j+1);
        if neighbor < subdomain,
            overlapIds = matchList{neighbor}(first:last) + int32( ones(numberShared,1) );
            overlap = forceVector{neighbor}(overlapIds);
            myNodalSum( sharedIds ) = myNodalSum( sharedIds ) + overlap;
        end
    end
    residual = myNodalSum - nodalSum{subdomain};
    assert( norm(residual) < 1.e-15 );
end
disp('Pass');
