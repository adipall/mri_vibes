function disp_gset=getdispr_g(cbmap,OTM,OutMap,fetimapdim,vr)
%function disp_gset=getdispr_g(cbmap,OTM,OutMap,fetimapdim,vr)
% computes the displacement in G space from the displacement
% in the reduced Craig-Bampton space.  Compare expandRmodel
dof_per_node = 8;
[vr_exp, nodes]=expandRmodel(cbmap,OTM,OutMap,vr);
disp_gset=zeros(fetimapdim,1);
for i=1:max(size(nodes))
  k=(nodes(i)-1)*dof_per_node;
  k1=(i-1)*6;       % 6?
  for cid=1:6
    disp_gset(k+cid)=vr_exp(k1+cid);
  end
end
