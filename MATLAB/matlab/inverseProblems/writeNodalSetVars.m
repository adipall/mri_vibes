load 2cubetest-out.mat
dispxdata=nvar01(nsnod01,:);
dispydata=nvar02(nsnod01,:);
dispzdata=nvar03(nsnod01,:);
pressuredata=nvar07(nsnod02,:);
new=exo_get('2cubetest.exo')        
nnset = size(new.Nodesets,2);

new.NodesetVars(1).NodesetVarData(1).Name = 'dispx';
new.NodesetVars(1).NodesetVarData(1).Data = dispxdata;
new.NodesetVars(1).NodesetVarData(2).Name = 'dispy';
new.NodesetVars(1).NodesetVarData(2).Data = dispydata;
new.NodesetVars(1).NodesetVarData(3).Name = 'dispz';
new.NodesetVars(1).NodesetVarData(3).Data = dispzdata;
new.NodesetVars(1).NodesetVarData(4).Name = 'APressure';
new.NodesetVars(1).NodesetVarData(4).Data = [];


new.NodesetVars(2).NodesetVarData(4).Name = 'APressure';
new.NodesetVars(2).NodesetVarData(4).Data = pressuredata;
%new.NodesetVars(2).NodesetVarData(2).Name = 'APressure';
new.NodesetVars(2).NodesetVarData(1).Data = [];
new.NodesetVars(2).NodesetVarData(2).Data = [];
%new.NodesetVars(2).NodesetVarData(3).Name = 'APressure';
new.NodesetVars(2).NodesetVarData(3).Data = [];
new.Time = [];
dt = 0.01;
for i=1:10000
  new.Time(i) = i*dt;
end
exo_put('2cubetest_new.exo',new)

