function UD = UO_Iwan(u,UD,params)

if isfield(UD,'init')
  [~, UD] = IMEX_Iwan(u,params) ;
  UD.init = 1;
else
  [~, UD] = IMEX_Iwan(u,UD) ;
end