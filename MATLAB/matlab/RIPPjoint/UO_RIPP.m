function UD = UO_RIPP(u,v,UD,RIPP)

if isfield(UD,'init')
  [~, UD] = RIPPjoint(u,v,RIPP) ;
  UD.init = 1;
else
  [~, UD] = RIPPjoint(u,v,UD) ;
end
