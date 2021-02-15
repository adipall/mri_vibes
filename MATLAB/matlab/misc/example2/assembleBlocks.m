% One path to determine dof[node_id] % from nodal_dof_per_block 
% leads through ...
% What happens if a node is in a pair of blocks with 
% different degrees of freedom? 
metaData,
eval('load olio');
act = 'aBlock=';
root = 'blk0';
number_blocks = nblks;  % each subdomain has the same number of blocks
clear nblks;
assert( 0<number_blocks && number_blocks < 10 );  % single digit
% todo: modify this to determine nodal_dof 
number_elements = zeros(number_blocks,1);
element_degree = zeros(number_blocks,1);
nodes = [];
number_nodes = 0;
for block = 1:number_blocks, % 1 <= number_blocks <= 9 
   eval( [act,root,int2str(block),';']);
   [m,n] = size(aBlock);
   number_elements(block) = m;
   element_degree(block) = n;
   p = m*n;
   if p > 0
       nodes( number_nodes+1: number_nodes + p ) = reshape( aBlock, p, 1);
   end
end
unique_nodes = unique(nodes);
end
