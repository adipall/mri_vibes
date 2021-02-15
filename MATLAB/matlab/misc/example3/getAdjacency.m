%function [matchListSize,matchlist,owner]= getAdjacency(numberSubdomains)
function [matchListSize,matchlist,owner]= getAdjacency(numberSubdomains)
matchListSize = cell(numberSubdomains,1);
matchlist = cell(numberSubdomains,1);
owner = cell(numberSubdomains,1);
for subdomain =1:numberSubdomains,
    [listSize, list, degree ] = getSubdomainAdjacency(subdomain);
    matchListSize{subdomain} = listSize;
    matchlist{subdomain} = list;
    owner{subdomain} = degree;
end
