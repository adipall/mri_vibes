function Q=t_window(x,N,lap)
% Q=t_window(x,N,lap)
% This function produces a data matrix from the data 'x'. Each column in
% the matrix is 'N' values long, overlapping 'lap' points, with the last
% one being zero-padded if need be.

% from T Edwards

nx       = length(x);
noverlap = lap;
nwind    = N;
ncol     = fix((nx-noverlap)/(nwind-noverlap));

colindex = 1 + (0:(ncol-1))*(nwind-noverlap);
rowindex = (1:nwind)';

% zero-padding data
if nx<(nwind+colindex(ncol)-1),
    x = 0; 
end

% output
Q    = rowindex(:,ones(1,ncol))+colindex(ones(nwind,1),:)-1;
Q    = zeros(nwind,ncol);
Q(:) = x(rowindex(:,ones(1,ncol))+colindex(ones(nwind,1),:)-1);
