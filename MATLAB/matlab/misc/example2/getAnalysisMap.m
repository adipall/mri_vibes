% function analysis_map = getAnalysisMap(subdomain)
function analysis_map = getAnalysisMap(subdomain)
analysis_map = [];
assert( 0<subdomain && subdomain < 10 ); % single digit
act = 'analysis_map = FetiMap_a_';
tag = int2str(subdomain-1);
suffix = '();'
eval([act,tag,,suffix]),
