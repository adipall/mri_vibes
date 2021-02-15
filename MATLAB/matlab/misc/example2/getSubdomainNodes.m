%function node_num_map = getSubdomainNodes(subdomain)
function node_num_map = getSubdomainNodes(subdomain)
node_num_map = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
act = 'load p';
tag = int2str(subdomain-1);
eval([act,tag]),
