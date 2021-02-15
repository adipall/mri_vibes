%function [forceVector,nodalSum,displacementVector]= getFeVectors(numberSubdomains)
function [forceVector,nodalSum,displacementVector]= getFeVectors(numberSubdomains)
forceVector = cell(numberSubdomains,1);
nodalSum = cell(numberSubdomains,1);
displacementVector = cell(numberSubdomains,1);
for subdomain =1:numberSubdomains,
    [force, summed, deform] = getSubdomainFeVectors(subdomain);
    forceVector{subdomain} = force;
    nodalSum{subdomain} = summed;
    displacementVector{subdomain} = deform;
end
