% function index =getTimeIndex(KF,time)
%                          KF=1   t=[0 2 4 6 ... (2*ncoeffs-2)]
%                          KF=2   t=[1 3 5 7 ... (2*ncoeffs-1)]
%                          KF=3   t=[2 4 6 8 ... (2*ncoeffs)]
%                          KF=4   t=[1 3 5 7 ... (2*ncoeffs-1)]
function index =getTimeIndex(KF,time)
if KF == 1      %    vt=2*ik;
    s=time/2;
    offset = 1;
elseif KF == 2  %    vt=2*ik+1;
    s=(time-1)/2;
    offset = 1;
elseif KF == 3  %    vt=2*ik;
    s =time/2;
    offset = 0;
 elseif KF == 4 %    vt=2*ik+1;
    s = (time-1)/2;
    offset = 1;
end 
index=fix(s)+offset;
