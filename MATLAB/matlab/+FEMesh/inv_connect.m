function [node_elem,nodes]=inv_connect(conn)
%
%
num_nodes_per_element=size(conn,2);
%
conn=conn';
[scon,idx]=sort(conn(:));

[nodes,idx2]=unique(scon,'last');
idx=idx-1;

idx=floor(idx/num_nodes_per_element)+1; % these are the elements corresponding to the nodes in scon
idx(~idx)=1;
%
node_elem=cell(length(nodes),1);
idxprev=1;
for i=1:length(nodes),
    node_elem{i}=idx(idxprev:idx2(i));
    idxprev=idx2(i)+1;
end