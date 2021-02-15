%function ik = getIdList(category,ncoeffs)
function ik = getIdList(category,ncoeffs)
if category == 1     %even-even
    ik=0:ncoeffs-1;
elseif category == 2 %even-odd
    ik=0:ncoeffs-1;
elseif category == 3 %odd-even
    ik=1:ncoeffs;
elseif category == 4 %odd-odd
    ik=0:ncoeffs-1;
end
ik = ik';
