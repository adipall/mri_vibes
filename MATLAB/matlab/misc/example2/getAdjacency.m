%function [matchListSize,matchlist,cardinality]= getAdjacency(numberSubdomains)
function [matchListSize,matchlist,cardinality]= getAdjacency(numberSubdomains)
matchListSize = cell(numberSubdomains,1);
matchlist = cell(numberSubdomains,1);
cardinality = cell(numberSubdomains,1);
for subdomain =1:numberSubdomains,
    [listSize, list, degree ] = getSubdomainAdjacency(subdomain);
    matchListSize{subdomain} = listSize;
    matchlist{subdomain} = list;
    cardinality{subdomain} = degree;
end
