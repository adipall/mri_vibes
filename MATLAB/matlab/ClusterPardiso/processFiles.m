procIDs = [0 1 2];
clear A
tic
tstart = tic;
numTerms = 0;
for i=1:length(procIDs)
  fileName = ['A_CP_' num2str(procIDs(i)) '.dat'];
  A{i} = rdCPfile(fileName);
  numTerms = numTerms + size(A{i},1);
end
telapsed = toc(tstart);
fprintf(1,'time to read files = %g\n', telapsed);
AA = zeros(numTerms, 3);
numTerms = 0;
for i=1:numProcs
  ii = numTerms+1:numTerms+size(A{i},1);
  AA(ii,:) = A{i};
  numTerms = numTerms + size(A{i},1);
end
clear A
A = spconvert(AA); % only upper triangule at this point
A = A + A' - diag(diag(A));
tic;
tstart = tic;
p = symamd(A);
R = chol(A(p,p));
telapsed = toc(tstart);
fprintf(1,'time for factorization = %g\n', telapsed);

