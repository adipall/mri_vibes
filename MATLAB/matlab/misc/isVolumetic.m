% function check = isVolumetic(fmap)
function check = isVolumetic(fmap)
p = size(fmap,1);
fmap = loadMaps(p);
eight = getDofPerNode();
five = eight - 3;
check = zeros(p,five);
for i=1:p,
    m = size(fmap{i},1);
    assert(rem(m,eight)== 0);
    for j=1:five,
       check(i,j)= max( fmap{i}(j+3:eight:m) );
    end
end
