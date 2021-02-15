% function nodes = getNodes(numberSubdomains)
function nodes = getNodes(numberSubdomains)
nodes = cell(numberSubdomains,1);
for subdomain =1:numberSubdomains,
     nodes{subdomain} = getSubdomainNodes(subdomain);
end
