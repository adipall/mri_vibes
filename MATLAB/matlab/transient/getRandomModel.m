%function [k,c,m] = getRandomModel()
function [k,c,m] = getRandomModel()
m = rand(); % small m <-> hard
k = rand();
c = rand();
while c*c > 4*m*k,
   c = c*.5;
end
