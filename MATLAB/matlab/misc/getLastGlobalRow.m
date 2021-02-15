% function last = getLastGlobalRow(globalRow)
function last = getLastGlobalRow(globalRow)
p = size(globalRow,1);
last = 0;
for i = 1:p,
  last=max([last,max(globalRow{i})]); 
end
last = last + 1;
