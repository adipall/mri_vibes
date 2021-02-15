%function [matchListSize,matchlist,cardinality]= getSubdomainAdjacency(subdomain)
numberSubdomains = int32(3);
[matchListSize,matchlist,cardinality]= getAdjacency(numberSubdomains);
for subdomain = 1:numberSubdomains,
    assert( 3 == size( matchListSize{subdomain},1 )   );
end
% subdomains node overlaps with its neighbors.
for subdomain = 1:numberSubdomains,
    assert( matchListSize{subdomain}(3) == size( matchlist{subdomain},1 )   );
end

for subdomain = 1:numberSubdomains,
    figure(subdomain);
    neighbor = 1;
    i = matchListSize{subdomain}(neighbor)+1;
    f = matchListSize{subdomain}(neighbor+1);
    plot( i:f,matchlist{subdomain}(i:f),'o'); hold on;
    i = f+1;
    neighbor = 1;
    f = matchListSize{subdomain}(3);
    plot( i:f,matchlist{subdomain}(i:f),'x'); hold off;
end
% zero-based
% values in the analysis set