%function [matchListSize,matchList,owner]= getSubdomainAdjacency(subdomain)
function [matchListSize,matchList,owner]= getSubdomainAdjacency(subdomain)
matchListSize = [];
matchList = [];
owner = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
tag = int2str(subdomain-1);
for species = 0:2, 
    switch species
        case 0
            root = 'matchListSize';
        case 1
            root = 'matchList';
        case 2
            root = 'owner';
        otherwise
            error('poach: illegal species');
    end % species
    underscore ='_';
    myFile = [root,underscore,tag];
    equals = '=';
    action = [root,equals,myFile,';'];
    eval( action );
    switch species
        case 0     
            matchListSize = int32( matchListSize ); 
        case 1
            matchList  = int32( matchList );
        case 2
            owner  = int32( owner );    
        otherwise
            error('illegal species');
    end % speciesj
end
