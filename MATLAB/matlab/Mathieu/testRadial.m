%%  test radial Mathieu functions
clear
nCoeffs = 5;
result = ones(1,4);
for category=1:4,
    [radii,q,value] = radialFunctionData(category);
    numberRadial = length(radii);
    numberQ = length(q);
    my = zeros(numberRadial,numberQ);
    for iq = 1:numberQ;
        [~,coefficients,~]=eig_Spm(category,q(iq),nCoeffs);
        for ku=1:numberRadial;
            nmax = 1;
            my(ku,iq)=Jpm(category,radii(ku),q(iq),coefficients,nmax);
        end
    end
    barmyNormalization = my(1,:)./value(1,:);
    err = value * diag(barmyNormalization) - my;
    result(category) = max(max(abs(err)));
end
if max(result) < 1.e-3
    disp('pass');
else
    disp('fail');
end
