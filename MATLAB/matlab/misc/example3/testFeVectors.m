% figure(1);
clear
numberSubdomains = 4;
subdAd = getSubdAd();
[forceVectorC, nodalSumC, displacementVectorC] = getFeVectors(subdomain);
for subdomain = 1:numberSubdomains,
    [forceVector, nodalSum, displacementVector] = getSubdomainFeVectors(subdomain);

    x = forceVectorC{subdomain} - forceVector;
    y = nodalSumC{subdomain} - nodalSum;
    z = displacementVectorC{subdomain} - displacementVector;
    assert(  norm( [x;y;z] ) == 0 );
end
disp('Passed');  % all assertions hold



