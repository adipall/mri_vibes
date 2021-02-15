% function coefficientIndices = getVt(catetory, numberCoefficients)
%  category=1   vt=[0 2 4 6 ... (2*ncoeffs-2)]
%  category=2   vt=[1 3 5 7 ... (2*ncoeffs-1)]
%  category=3   vt=[2 4 6 8 ... (2*ncoeffs)]
%  category=4   vt=[1 3 5 7 ... (2*ncoeffs-1)]
function vt = getVt(category,numberCoefficients)
ik = getIdList(category,numberCoefficients);
if category == 1 || category == 3, %even-even v odd-even 
    vt=2*ik;
else      % category == 2v4  even-odd v odd-odd
    vt=2*ik+1;
end
