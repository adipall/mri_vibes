% function value=extract_one_value(KF,t,vncoeffs)
function value=extract_one_value(KF,t,vncoeffs)
index = getTimeIndex(KF,t);
value=vncoeffs(index);
