clear
%% test series
parameter = .1;              % q
numberModes = 5;             % n
numberTerms = numberModes+3; % m
numberCoordinates = 51;
angles = linspace(0,pi/2,numberCoordinates)';
testResults = zeros(4,1);

for category=1:4
    [~,coefficients,indices]=eig_Spm(category,parameter,numberTerms);
    if category == 1 || category == 2,
        trigonometry = 'cosine';
    else
        trigonometry = 'sine';
    end
    my=series(angles,coefficients,indices,trigonometry);
    values = zeros(numberTerms,numberCoordinates);
    for i=1:numberCoordinates,
        values(:,i)=Spm(category,angles(i),coefficients,numberTerms);
    end
    testResults(category)=norm(my - values');
end
if norm(testResults,inf) < 1.e-12,
    disp('pass');
else
    disp('fail');
end

