function A = rdFile(fileName)
%function A = rdFile(fileName)
% reads Cluster Pardiso text file containing matrix information
fileID = fopen(fileName, 'r');
range = cell2mat(textscan(fileID, '%f', 2, 'Delimiter', '\n'));
firstRow = range(1);
lastRow = range(2);
numRows = lastRow - firstRow + 1;
rowBegin = cell2mat(textscan(fileID, '%f', numRows+1, 'Delimiter', '\n'));
numTerms = rowBegin(end) - 1;
numCols = 2;
B = cell2mat(textscan(fileID, '%f %f', numTerms, 'Delimiter', '\n'));
A = zeros(numTerms, 3);
for i=1:numRows
  row = firstRow + i - 1;
  for j=rowBegin(i):rowBegin(i+1)-1
    col = B(j,1);
    if (j > rowBegin(i))
      if (col <= B(j-1,1))
        error('columns not in ascending order');
      end
    end
    A(j,:) = [row col B(j,2)];
  end
end
fclose(fileID);
