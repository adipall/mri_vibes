numberSubdomains = 3;
nodes = getNodes(numberSubdomains);
allConnected = numberSubdomains*( numberSubdomains + 1 )/2;
pairs = zeros(allConnected,2);
overlap = cell(numberSubdomains,1);
i = 0;
for subdomain = 1:numberSubdomains-1,
    for adjacentSubdomain  = subdomain+1:numberSubdomains,
        i = i + 1;
        pairs(i,:) = [subdomain, adjacentSubdomain];
        overlap{i} = intersect( nodes{subdomain}, nodes{adjacentSubdomain})';
    end
end
