%function [forceVector, nodalSum, displacementVector] = getSubdomainFeVectors(subdomain)
function [forceVector, nodalSum, displacementVector] = getSubdomainFeVectors(subdomain)
forceVector = [];
nodalSum = [];
displacementVector = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
vectorNames = getVectorNames();
numSpecies = size(vectorNames,1);
tag = int2str(subdomain-1);
last = numSpecies-1;
for species = 0:last,
    root = vectorNames{species+1};
    myFile = strcat( root, tag );
    run(myFile);
    switch species
        case 0
           % forceVector
        case 1
            % nodalSum
        case 2             
%             displacementVector    
        otherwise
            error('illegal species');
    end % speciesj
end
