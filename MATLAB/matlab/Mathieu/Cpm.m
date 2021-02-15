%   function correlationFactors = Cpm(category,mv,mvv,nmax)
%   mv = matrix of expansion coefficients from 'eig_Spm'
%   mvv= matrix of expansion coefficients from 'eig_Spm'
%   mv and mvv are square matrices
%   nmax= maximum order
%   category:(even even,even odd,odd even,odd odd)
function correlationFactors = Cpm(category,mv,mvv,nmax)
if  size(mv,2) ~= size(mvv,2)
    correlationFactors =[];
    return;
end
nCoeffs = size(mv,2);
if  (category < 1) || (category > 4) || (nmax <= 0) || (nmax > nCoeffs)
    correlationFactors =[];
    return;
end
correlationFactors=zeros(1,nmax);
if category == 1, % even-even
    for k=1:nmax,
        Apm=mv(:,k);
        AApm=mvv(:,k);
        correlationFactors(k)=...
            2*Apm(1)*AApm(1)+sum(Apm(2:nCoeffs).*AApm(2:nCoeffs));
    end
else % even odd,odd even,odd odd
    for k=1:nmax,
        correlationFactors(k)= sum(mv(:,k).*mvv(:,k));
    end
end
correlationFactors = correlationFactors*pi;
