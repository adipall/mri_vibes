function [dispgr,nodes]=expandRmodel( cbmap, OTM, OutMap, vr )
%function [dispgr,nodes]=expandRmodel( cbmap, OTM, OutMap, vr )
%
% expands a vector in the reduced, craig-bampton space into the
% full physical space.
% cbmap - map to interface dofs. Output into cbr.m
% OTM - Output transfer matrix. also in cbr.m
% OutMap - map to interior (and perhaps interface) nodes in output.
% vr - the reduced space vector.
%     vr(1:numeig) is the amplitude of the fixed interface modes
%     vr(numeig:end) is the amplitude of the constraint modes (physical
%                        degrees of freedom).
% results are output sorted by node number. 6 dofs per node are output.

% $Revision$
% $Date$

nodes=unique([cbmap(:,1)' OutMap]);  % sorted list of nodes output
nr=max(size(vr));                    % number dofs in vr
nc=size(cbmap,1);                    % number of constraint dofs
nmodes=nr-nc;                        % number of fixed interface modes
dispgr=zeros(max(size(nodes))*6,1);
ur=OTM*vr;                           % compute vector on OTM space, ur

% store components from OTM space.
for i=1:size(OutMap,2)
  k=find(nodes==OutMap(i));
  for cid=1:6
    k2=(k-1)*6+cid;
    k1=(i-1)*6+cid;
    dispgr(k2)=ur(k1);
  end
end

% transfer interface dofs directly
for i=1:nc
  n=cbmap(i,1);
  cid=cbmap(i,2);
  k=find(nodes==n);
  k2=(k-1)*6+cid;
  dispgr(k2)=vr(i+nmodes);
end
