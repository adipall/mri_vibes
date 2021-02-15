%function gid = getSubdomainNodes2(subdomain)
function gid = getSubdomainNodes2(subdomain)
gid = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
tag = int2str(subdomain-1);
act = 'olio_gid_';
eval([act,tag]),
doublePrecision = gid;
gid = int32( doublePrecision );

