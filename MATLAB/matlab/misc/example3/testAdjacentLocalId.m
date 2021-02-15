clear;
numberSubdomains = 4;
subdAd = getSubdAd();
adjacentLocalId = getAdjacentLocalId();
assert( size( subdAd, 1) == numberSubdomains );
assert( size( adjacentLocalId, 1) == numberSubdomains );
for subdomain = 1:numberSubdomains,
    assert( size( subdAd{subdomain}, 1) == size( adjacentLocalId{subdomain}, 1))
end
for subdomain = 1:numberSubdomains,
    n = size( subdAd{subdomain}, 1);
    for i=1:n,
        neighbor = subdAd{subdomain}(i);
        j = adjacentLocalId{subdomain}(i);
        k = subdAd{neighbor}(j);
        assert( k == subdomain );
    end
end

[matchListSize,matchList,owner]= getAdjacency(numberSubdomains);
for subdomain = 1:numberSubdomains,
    numberNeighbors = size(matchListSize{subdomain},1) - 1;
    for i = 1:numberNeighbors,
        neighbor = subdAd{subdomain}(i);
        start = matchListSize{subdomain}(i) + 1;
        stop = matchListSize{subdomain}(i+1);
        left = stop - start + 1;
        j = adjacentLocalId{subdomain}(i);
        hi = matchListSize{neighbor}(j) + 1;
        ho = matchListSize{neighbor}(j+1);
        right = ho - hi + 1;
        assert( left == right );
    end
end