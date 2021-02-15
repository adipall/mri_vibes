% figure(1);
clear
numberSubdomains = 4;
subdAd = getSubdAd();
[matchListSize,matchList,owner]= getAdjacency(numberSubdomains);
dots = zeros(numberSubdomains,1);
for subdomain = 1:numberSubdomains,
    [forceVector, nodalSum, displacementVector] = getSubdomainFeVectors(subdomain);
    dots(subdomain) = nodalSum'*displacementVector;
end