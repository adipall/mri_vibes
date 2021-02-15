%%  chen 1994  ... I don't know either.
clear
nCoeffs = 30; % 20
result = ones(1,4);
numberRadial = 40; % 100
Ro = 2;           // R is 2
dead = 0;
notAvailable = -1;
radii = linspace(0,Ro,numberRadial);
perm54 = [1,3,2,4];
for category=3:3,
    row = perm54(category);
   row = 5;
    k = [];
    for column = 1:7,
        ko = kTableRis2(row,column);
        if ko ~= dead,
            k = [k; ko];
        end
    end
    q = (.5*k).^2;
    numberQ = length(q);
    my = zeros(numberRadial,numberQ);
    for iq = 1:numberQ;
        [~,coefficients,~]=eig_Spm(category,q(iq),nCoeffs);
        for ku=1:numberRadial;
            nmax = 1;
            my(ku,iq)=Jpm(category,radii(ku),q(iq),coefficients,nmax);
        end
    end
    plot(my),
end
