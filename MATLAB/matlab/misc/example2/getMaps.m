% function maps = getMaps(numberSubdomains)
function maps = getMaps(numberSubdomains)
maps = cell(numberSubdomains,1);
for subdomain =1:numberSubdomains,
     maps{subdomain} = getSubdomainMaps(subdomain);
end
