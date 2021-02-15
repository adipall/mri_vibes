function disp_gset=getdispg(FetiMap_a,dispAset)
%function disp_gset=getdispg(FetiMap_a,dispAset)
% computes the displacement in G space from the displacement
% in A-set space.
ng=size(FetiMap_a,1);
disp_gset=zeros(ng,1);

for i=1:ng
  k=FetiMap_a(i);
  if (k>0)
    disp_gset(i)=dispAset(k);
  end
end
