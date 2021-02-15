%function [matchListSize,matchlist,cardinality]= getSubdomainAdjacency(subdomain)
function [matchListSize,matchlist,cardinality]= getSubdomainAdjacency(subdomain)
matchListSize = [];
matchlist = [];
cardinality = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
tag = int2str(subdomain-1);
for species = 0:2, 
    switch species
        case 0
            root = 'matchListSize';
        case 1
            root = 'matchlist';
        case 2
            root = 'cardinality';
        otherwise
            error('poach: illegal species');
    end % species
    myFile = [root,tag];
    switch species
        case 0
            matchListSize = load(myFile);
            matchListSize = int32( matchListSize(:,2) ); 
        case 1
            matchlist  = load(myFile);
            matchlist  = int32( matchlist(:,2) );
        case 2
            cardinality  = load(myFile);
            cardinality  = int32( cardinality(:,2) );
        otherwise
            error('illegal species');
    end % speciesj
end
