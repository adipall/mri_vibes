% function coords = quad(rMin, rMax, thetaMin, thetaMax)
function coords = quad(rMin, rMax, thetaMin, thetaMax)
index = 1;
coords(index,1) = rMin * cos(thetaMin);
coords(index,2) = rMin * sin(thetaMin);
index = index + 1;

coords(index,1) = rMin * cos(thetaMax);
coords(index,2) = rMin * sin(thetaMax);
index = index + 1;

coords(index,1) = rMax * cos(thetaMax);
coords(index,2) = rMax * sin(thetaMax);
index = index + 1;

coords(index,1) = rMax * cos(thetaMin);
coords(index,2) = rMax * sin(thetaMin);
index = index + 1;

coords(index,1) = rMin * cos(thetaMin);
coords(index,2) = rMin * sin(thetaMin);
