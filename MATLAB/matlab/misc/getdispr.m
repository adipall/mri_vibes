function dispr=getdispr(cbmap,disp_gset)
%function dispr=getdispr(cbmap,disp_gset)
%
% computes the reduced space vector (u_{interface}) given
% the displacements in global physical space (G-set).
dof_per_node = getDofPerNode();
k=(cbmap(:,1)-1)*dof_per_node + cbmap(:,2);
dispr=disp_gset(k);
