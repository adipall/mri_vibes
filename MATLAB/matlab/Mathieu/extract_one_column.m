% function vec=extract_one_column(KF,t,mvncoeffs)
function vec=extract_one_column(KF,t,mvncoeffs)
index = getTimeIndex(KF,t);
vec=mvncoeffs(:,index);
