%function matrix=getMatrix(category,parameter,numberCoefficients)
function matrix=getMatrix(category,parameter,numberCoefficients)
matrix = [];
if category < 1 || category  > 4 ||  parameter <= 0 || numberCoefficients < 3
   return;
end
ik = getIdList(category,numberCoefficients);
upDiagonal = parameter*ones(1,numberCoefficients-1); 
if category == 1
    lowDiagonal = parameter*[2, ones(1,numberCoefficients-2)];
    diagonal = (2*ik).^2; 
elseif category == 2  
    lowDiagonal = upDiagonal;
    diagonal = (2*ik+1).^2;
    diagonal(1) = diagonal(1) + parameter;
elseif category == 3  
    lowDiagonal = upDiagonal;
    diagonal = (2*ik).^2;
 elseif category == 4  
    lowDiagonal = upDiagonal;
    diagonal = (2*ik+1).^2;
    diagonal(1) = diagonal(1) - parameter;
end 
matrix = diag(diagonal) + diag( lowDiagonal,-1) + diag(upDiagonal, 1);
