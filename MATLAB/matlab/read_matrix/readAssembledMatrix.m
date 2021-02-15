function A = readAssembledMatrix(numProc)
%function A = readAssembledMatrix(numProc)
% numProc = number of processors for assembled matrix
numTermsAll = 0;
for i=1:numProc
  mat = ['assembledMatrix' num2str(i-1)];
  fname = [mat '.dat;'];
  a = ['load ' fname];
  eval(a);
  b = ['numTerms = size(' mat ',1);'];
  eval(b)
  numTermsAll = numTermsAll + numTerms;
end
AA = zeros(numTermsAll,3);
numTermsAll = 0;
for i=1:numProc
  mat = ['assembledMatrix' num2str(i-1)];
  b = ['numTerms = size(' mat ',1);'];
  eval(b)
  ii = numTermsAll+1:numTermsAll+numTerms;
  AA(ii,:) = eval(mat);
  numTermsAll = numTermsAll + numTerms;  
end
A = spconvert(AA);
