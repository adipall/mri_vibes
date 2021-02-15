function [domains,numberSharedNodes,sharedGids,currentLids,adjacentLids]= getSharedNodes(gids)
%function [domains,numberSharedNodes,sharedGids,currentLids,adjacentLids]= getSharedNodes(gids)
% example 
%  gids = loadGids(root,p);
%  [d,n,s,c,a] = getSharedNodes( gids );
%    
% See also setNodeMaps, overlaps
assert( size( gids, 2 ) == 1 );
p = size( gids , 1 );
sharedGids = cell(p,1);
currentLids = cell(p,1);
adjacentLids = cell(p,1);
domains = cell(p,1);
numberSharedNodes = cell(p,1);
for subdomain = 1:p,
    m = 0;
    neighbor = 0;
    for adjDomain = 1:p,
    if( adjDomain ~= subdomain )
        [shared, current, adjacent] = intersect(gids{subdomain},gids{adjDomain});
        if ~isempty(shared)
            neighbor = neighbor + 1;
            domains{subdomain}(neighbor) = adjDomain;
            k = size(shared,1);
            numberSharedNodes{subdomain}(neighbor) = k;
            sharedGids{subdomain}(m+1:m+k) = shared;
            currentLids{subdomain}(m+1:m+k) = current;
            adjacentLids{subdomain}(m+1:m+k) = adjacent;
        end
    end
    end
end
